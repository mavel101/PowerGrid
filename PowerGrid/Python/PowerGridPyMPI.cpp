
/*
(C) Copyright 2015-2016 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/*****************************************************************************

    File Name   [PowerGridSenseMPI.cpp]

    Synopsis    [PowerGrid reconstruction executable supporting ISMRMRD format
                                as input.]

    Description [This reconstruction supports 2D and 3D reconstructions with
                                time segmentation for field correction. This
 support is experimental.]

    Revision    [0.2.0; Alex Cerjanic, BIOE UIUC]

    Date        [2019/03/31]

 *****************************************************************************/

// //Project headers.
#include "PowerGrid.h"
#include "processIsmrmrd.hpp"
#include "processNIFTI.hpp"
#include <boost/program_options.hpp>
#include <chrono>

#include "pybind11/pybind11.h"
#include "pybind11/numpy.h"
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include <pybind11/functional.h>
#include <pybind11/chrono.h>

namespace py = pybind11;
namespace bmpi = boost::mpi;

typedef std::vector<std::complex<float>> cmplx_vec;
typedef std::vector<std::vector<std::complex<float>>> cmplx_vec2d;

py::dict PowerGridSenseMPI(std::string inFile, std::string outFile, int nx, int ny, int nz, int nShots, std::string  TSInterp,
                     std::string FourierTrans, int timesegs, bool ts_adapt, double beta, int niter, int regDims) {

  // save image data and metadata in a dict
  py::dict imgs;

  uword Nx = nx;
  uword Ny = ny;
  uword Nz = nz;
  uword NShots = nShots;
  uword L = timesegs;
  uword NIter = niter;
  uword dims2penalize = regDims;
  uword FtType = 1, type = 1;

  bmpi::environment env(true);
  bmpi::communicator world;

    if (FourierTrans.compare("DFT") == 0) {
        FtType = 2;
    
    } else if (FourierTrans.compare("NUFFT") == 0) {
      
      FtType = 1;
      if (TSInterp.compare("hanning") == 0) {
        type = 1;
      } else if (TSInterp.compare("minmax") == 0) {
        type = 2;
      } else if (TSInterp.compare("histo") == 0) {
        type = 3;
      } else {
        std::cout << "Did not recognize temporal interpolator selection. " << std::endl
                  << "Acceptable values are hanning or minmax."            << std::endl;
        return imgs;
      }
    } else if (FourierTrans.compare("DFTGrads") == 0) {
      FtType = 3;
    } else {
      std::cout << "Did not recognize Fourier transform selection. " << std::endl
                << "Acceptable values are DFT or NUFFT."             << std::endl;
      return imgs;
    }


  ISMRMRD::Dataset *d;
  ISMRMRD::IsmrmrdHeader hdr;
	acqTracking *acqTrack;
  Col<float> FM;
  Col<std::complex<float>> sen;
 	processISMRMRDInput<float>(inFile, d, hdr, FM, sen, acqTrack);

  uword numAcq = d->getNumberOfAcquisitions();

  // Grab first acquisition to get parameters (We assume all subsequent
  // acquisitions will be similar).
  ISMRMRD::Acquisition acq;
  d->readAcquisition(numAcq-1, acq); // take last acquisition as first acquisition is usually a noise scan - mv
  uword nro = acq.number_of_samples();
  uword nc = acq.active_channels();

	// Handle Nx, Ny, Nz
	if(Nx==0) {
		Nx = hdr.encoding[0].encodedSpace.matrixSize.x;
	} 
	if(Ny==0) {
		Ny = hdr.encoding[0].encodedSpace.matrixSize.y;
	} 
	if(Nz==0) {
		Nz = hdr.encoding[0].encodedSpace.matrixSize.z;
	} 

  Col<float> ix, iy, iz;
  initImageSpaceCoords(ix, iy, iz, Nx, Ny, Nz);


  // Check and abort if we have more than one encoding space (No Navigators for
  // now).
  std::cout << "hdr.encoding.size() = " << hdr.encoding.size() << std::endl;
  if (hdr.encoding.size() != 1) {
    std::cout << "There are " << hdr.encoding.size()
              << " encoding spaces in this file" << std::endl;
    std::cout
        << "This recon does not handle more than one encoding space. Aborting."
        << std::endl;
    return imgs;
  }

  int NShotMax  = hdr.encoding[0].encodingLimits.kspace_encoding_step_1->maximum + 1;
  int NParMax   = hdr.encoding[0].encodingLimits.kspace_encoding_step_2->maximum + 1;
  int NSliceMax = hdr.encoding[0].encodingLimits.slice->maximum + 1;
  int NSetMax   = hdr.encoding[0].encodingLimits.set->maximum + 1;
  int NRepMax   = hdr.encoding[0].encodingLimits.repetition->maximum + 1;
  int NAvgMax   = hdr.encoding[0].encodingLimits.average->maximum + 1;
  int NSegMax   = hdr.encoding[0].encodingLimits.segment->maximum + 1;
  int NEchoMax  = hdr.encoding[0].encodingLimits.contrast->maximum + 1;
  int NPhaseMax = hdr.encoding[0].encodingLimits.phase->maximum + 1;

  std::cout << "NParMax = "   << NParMax << std::endl;
  std::cout << "NShotMax = "  << NShotMax << std::endl;
  std::cout << "NSliceMax = " << NSliceMax << std::endl;
  std::cout << "NSetMax = "   << NSetMax << std::endl;
  std::cout << "NRepMax = "   << NRepMax << std::endl;

  std::cout << "NAvgMax = "   << NAvgMax << std::endl;

  std::cout << "NEchoMax = "  << NEchoMax << std::endl;

  std::cout << "NPhaseMax = " << NPhaseMax << std::endl;

  std::cout << "NSegMax = "   << NSegMax << std::endl;


  std::cout << "About to loop through the counters and the file"
            << std::endl;

  py::list shapes;
  shapes.append(NSliceMax);
  shapes.append(NPhaseMax);
  shapes.append(NEchoMax);
  shapes.append(NAvgMax);
  shapes.append(NRepMax);
  shapes.append(Nz);
  shapes.append(Ny);
  shapes.append(Nx);
  imgs["shapes"] = shapes;

  // std::string baseFilename = "img";
  std::string filename;
  // if (!outputImageFilePath.empty() && *outputImageFilePath.rbegin() != '/') {
  //   	outputImageFilePath += '/';
	// }

  // Figure out what slices to work on in this rank
  Col<uword> sliceList;
  Col<uword> phaseList;
  Col<uword> echoList;
  Col<uword> avgList;
  Col<uword> repList;

  uword numTasks = NSliceMax * NPhaseMax * NEchoMax * NAvgMax * NRepMax;

  sliceList.resize(numTasks);
  phaseList.resize(numTasks);
  echoList.resize(numTasks);
  avgList.resize(numTasks);
  repList.resize(numTasks);

  // Now we make lists of the tasks such that we are mapping out all the tasks
  for (uword ii = 0; ii < NSliceMax; ii++) {   // slice loop
    for (uword jj = 0; jj < NPhaseMax; jj++) { // phase Loop
      for (uword kk = 0; kk < NEchoMax; kk++) {  // echo loop
        for (uword ll = 0; ll < NAvgMax; ll++) { // avg Loop
          for (uword mm = 0; mm < NRepMax; mm++) {   // rep loop
            sliceList( ii * NPhaseMax * NEchoMax * NAvgMax * NRepMax + jj * NEchoMax * NAvgMax * NRepMax + kk * NAvgMax * NRepMax + ll * NRepMax + mm) = ii;
            phaseList( ii * NPhaseMax * NEchoMax * NAvgMax * NRepMax + jj * NEchoMax * NAvgMax * NRepMax + kk * NAvgMax * NRepMax + ll * NRepMax + mm) = jj;
            echoList ( ii * NPhaseMax * NEchoMax * NAvgMax * NRepMax + jj * NEchoMax * NAvgMax * NRepMax + kk * NAvgMax * NRepMax + ll * NRepMax + mm) = kk;
            avgList  ( ii * NPhaseMax * NEchoMax * NAvgMax * NRepMax + jj * NEchoMax * NAvgMax * NRepMax + kk * NAvgMax * NRepMax + ll * NRepMax + mm) = ll; 
            repList  ( ii * NPhaseMax * NEchoMax * NAvgMax * NRepMax + jj * NEchoMax * NAvgMax * NRepMax + kk * NAvgMax * NRepMax + ll * NRepMax + mm) = mm;  
          }
        }
      }
    }
  }

  // Start by creating a 2D vector (vector of vectors) to hold the MPI rank ->
  // task mapping
  vector<uword> init(0);
  std::vector<std::vector<uword> > *taskList = new std::vector<std::vector<uword>>(world.size(), init);
  uword process = 0;
  for (uword task = 0; task < numTasks; task++) {
    (*taskList)[process].push_back(task);

    process++; 
    if (process == world.size()) {
      process = 0;
    }
  }
  uword taskIndex;
  std::cout << " Number of tasks = " << numTasks << std::endl;

  #ifdef OPENACC_GPU
  // Spread out our tasks across GPUs.
  uword gpunum;
  uword ngpus = acc_get_num_devices( acc_device_nvidia );
  if (ngpus) {
    gpunum = world.rank() % ngpus;
    acc_set_device_num( gpunum, acc_device_nvidia );
  } else {
    acc_set_device_type( acc_device_host );
  }
  #endif

  std::cout << "Rank = " << world.rank() << std::endl;

  uword NSlice, NRep, NAvg, NEcho, NPhase;
  uword L_save;
  double FM_range;
  double FM_range_ref;

  // define 2D vector for image data
  int n_imgs = (*taskList)[world.rank()].size();
  cmplx_vec2d img_data2d(n_imgs, cmplx_vec (Nx*Ny*Nz, 0)); // create 2D vector to save images thread safety

  for (uword ii = 0; ii < (*taskList)[world.rank()].size(); ii++) {

    // check for keyboard interrupt
    if (PyErr_CheckSignals() != 0)
      throw py::error_already_set();

    taskIndex = (*taskList)[world.rank()].at(ii);

    NSlice = sliceList(taskIndex);
    NRep   = repList(taskIndex);
    NAvg   = avgList(taskIndex);
    NEcho  = echoList(taskIndex);
    NPhase = phaseList(taskIndex);

    std::cout << "NSlice = " << NSlice << std::endl;
    std::cout << "NRep   = " << NRep   << std::endl;
    std::cout << "NAvg   = " << NAvg   << std::endl;
    std::cout << "NEcho  = " << NEcho  << std::endl;
    std::cout << "NPhase = " << NPhase << std::endl;
    

    //Col<float> FM;
    Col<float> fmSlice;
    //Col<std::complex<float>> sen;
 	  Col<std::complex<float>> senSlice;
    Col<float> kx(nro), ky(nro), kz(nro), tvec(nro);
    Col<std::complex<float>> data(nro * nc);
    Col<std::complex<float>> ImageTemp(Nx * Ny * Nz);

                      filename = outFile + "_" + "Slice" + std::to_string(NSlice) +
          								"_" + "Rep" + std::to_string(NRep) + "_" + "Avg" + std::to_string(NAvg) +
          								"_" + "Echo" + std::to_string(NEcho) + "_" + "Phase" + std::to_string(NPhase);
 
	                      senSlice = getISMRMRDCompleteSENSEMap<std::complex<float>>(d, sen, NSlice, Nx*Ny*Nz);
						            fmSlice = getISMRMRDCompleteFieldMap<float>(d, FM, NSlice, (uword) (Nx*Ny*Nz));
	                      getCompleteISMRMRDAcqData<float>(d, acqTrack, NSlice, NRep, NAvg, NEcho, NPhase, data, kx, ky,
			                    kz, tvec);

                      // Adapt number of time segments based on the range of the field map
                      if (ts_adapt){
                        L_save = L;
                        if (NSlice==0){
                          FM_range_ref = arma::max(arma::vectorise(fmSlice)) - arma::min(arma::vectorise(fmSlice));
                        }
                        else{
                          FM_range = arma::max(arma::vectorise(fmSlice)) - arma::min(arma::vectorise(fmSlice));
                          L = (int) (L*sqrt(FM_range/FM_range_ref));
                          std::cout << "Adapting time segments to L = " << L << " based on Field Map range." << std::endl; 
                        }
                      }

	                    std::cout << "Number of elements in kx = " << kx.n_rows << std::endl;
	                    std::cout << "Number of elements in ky = " << ky.n_rows << std::endl;
	                    std::cout << "Number of elements in kz = " << kz.n_rows << std::endl;
	                    std::cout << "Number of rows in data = " << data.n_rows << std::endl;
	                    std::cout << "Number of columns in data = " << data.n_cols << std::endl;

	                    QuadPenalty<float> R(Nx, Ny, Nz, beta, dims2penalize);

                      if (FtType == 1) {
	                      Gnufft<float> G(kx.n_rows, (float) 2.0, Nx, Ny, Nz, kx, ky, kz, ix,
			                    iy, iz);
	                      TimeSegmentation<float, Gnufft<float>> A(G, fmSlice, tvec, kx.n_rows, Nx*Ny*Nz, L, type, NShots);
                        SENSE<float, TimeSegmentation<float, Gnufft<float>>> Sg(A, senSlice, kx.n_rows, Nx*Ny*Nz, nc);
	                      ImageTemp = reconSolve<float, SENSE<float, TimeSegmentation<float, Gnufft<float>>>,
			                    QuadPenalty<float>>(data, Sg, R, kx, ky, kz, Nx,
			                    Ny, Nz, tvec, NIter);
                      } else if (FtType == 2) {
                        Gdft<float> A(kx.n_rows, Nx*Ny*Nz,kx,ky,kz,ix,iy,iz,fmSlice,tvec);
	                      SENSE<float, Gdft<float>> Sg(A, senSlice, kx.n_rows, Nx*Ny*Nz, nc);
	                      ImageTemp = reconSolve<float, SENSE<float, Gdft<float>>,
	                            QuadPenalty<float>>(data, Sg, R, kx, ky, kz, Nx,
                              Ny, Nz, tvec, NIter);
                      } else if (FtType == 3) {
                        GdftR2<float> A(kx.n_rows, Nx*Ny*Nz,kx,ky,kz,ix,iy,iz,fmSlice,tvec,Nx,Ny,Nz);
	                      SENSE<float, GdftR2<float>> Sg(A, senSlice, kx.n_rows, Nx*Ny*Nz, nc);
	                      ImageTemp = reconSolve<float, SENSE<float, GdftR2<float>>,
	                            QuadPenalty<float>>(data, Sg, R, kx, ky, kz, Nx,
                              Ny, Nz, tvec, NIter);
                      }

                    if (!outFile.empty())
                      writeNiftiMagPhsImage<float>(filename,ImageTemp,Nx,Ny,Nz);
                      
                    // save data for Python conversion
                    for(int j = 0; j < Nx * Ny * Nz; j++)
                      img_data2d[ii][j] = static_cast<std::complex<float>>(ImageTemp(j));

                    // set L back to original value
                    L = L_save;

                    // check for keyboard interrupt
                    if (PyErr_CheckSignals() != 0)
                      throw py::error_already_set();

                    
    }

  // copy image data to pybind dict
  imgs["img_data"] = img_data2d;

  // Close ISMRMRD::Dataset, hdr, and acqTrack
	closeISMRMRDData(d,hdr,acqTrack);

  return imgs;
}


PYBIND11_MODULE(PowerGridPyMPI, m) {

  // Expose the function PowerGridSenseMPI().
  m.def("PowerGridSenseMPI", &PowerGridSenseMPI, 
              "Executes PowerGridSenseMPI with specified parameters.\n\\
                Parameters\n\\
                    inFile : string\n\\
                        input ISMRMRD Raw Data file (incl path)\n\\
                    outFile : string\n\\
                        output path for Nifti output file (incl path) - default not saving Images into Nifti file\n\\
                    FourierTrans : string\n\\
                        Fourier Transform (NUFFT (default), DFT or DFTGrads)\n\\
                    timesegs : int\n\\
                        Number of time segments for B0 correction - 0: no B0 correction, -1 (default): segments will be determined from readout duration\n\\
                    niter: int\n\\
                        Number of CG iterations (default=10)\n\\
                    TSInterp: string\n\\
                        Field Correction Interpolator (histo (default), hanning or minmax)\n\\
                    beta: float\n\\
                        Spatial regularization penalty weight (default=0)\n\\
                    regDims: int\n\\
                        Dimensions to apply regularization to (2 or 3 (default))\n\\
                    nShots: int\n\\
                        Number of Shots (default: 1)\n\\
                    nx, ny, nz : int\n\\
                        Image size in x,y,z. Default is 0, then sizes are read from ISRMRD header encoded space.\n\\
                Returns\n\\
                    Dict containing the image vector and corresponding shapes. Image shape can be regained doing:\n\\
                    np.asarray(dict[\"img_data\"]).reshape(dict[\"shapes\"])\n",
                py::arg("inFile"), py::arg("outFile")="", py::arg("nx")=0, py::arg("ny")=0, py::arg("nz")=0, py::arg("nShots")=1, py::arg("TSInterp")="histo",
                 py::arg("FourierTrans")="NUFFT", py::arg("timesegs")=-1, py::arg("ts_adapt")=false, py::arg("beta")=0.0, py::arg("niter")=10, py::arg("regDims")=3
        );
}