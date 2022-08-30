
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
#include "cnpy.h"

namespace po = boost::program_options;
namespace bmpi = boost::mpi;

typedef std::vector<std::complex<float>> cmplx_vec;

int main(int argc, char **argv) {
  std::string rawDataFilePath, outputImageFilePath, senseMapFilePath,
      fieldMapFilePath, precisionString,
      rawDataNavFilePath;

  bmpi::environment env(argc, argv, true);
  bmpi::communicator world;

  uword Nx, Ny, Nz, NIter = 10;
  //uword ;
  double beta = 0.0;
  uword dims2penalize = 3;
  bool writeNifti;
  po::options_description desc("Allowed options");
  desc.add_options()("help,h", "produce help message")(
      "inputData,i", po::value<std::string>(&rawDataFilePath)->required(),
      "input ISMRMRD Raw Data file")
 			("outputImage,o", po::value<std::string>(&outputImageFilePath)->required(), "output file path for Numpy files and NIFTIimages")
			("writeNifti,w", po::bool_switch(&writeNifti)->default_value(false), "Write each image in own Nifti file if selected.")
      ("Nx,x", po::value<uword>(&Nx), "Image size in X")
			("Ny,y", po::value<uword>(&Ny), "Image size in Y")
			("Nz,z", po::value<uword>(&Nz), "Image size in Z")
          ("Beta,B", po::value<double>(&beta), "Spatial regularization penalty weight")
          ("CGIterations,n", po::value<uword>(&NIter), "Number of preconditioned conjugate gradient interations for main solver")
          ("Dims2Penalize,D", po::value<uword>(&dims2penalize), "Dimensions to apply regularization to (2 or 3).");


  po::variables_map vm;

  try {

    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help")) {
      std::cout << desc << std::endl;
      return 1;
    }

  } catch (boost::program_options::error &e) {
    std::cerr << "Error: " << e.what() << std::endl;
    std::cout << desc << std::endl;
    return 1;
  }


  ISMRMRD::Dataset *d;
  ISMRMRD::IsmrmrdHeader hdr;
	acqTracking *acqTrack;
  Col<float> FM;
  Col<std::complex<float>> sen;
 	processISMRMRDInput<float>(rawDataFilePath, d, hdr, FM, sen, acqTrack);

  uword numAcq = d->getNumberOfAcquisitions();

  // Grab first acquisition to get parameters (We assume all subsequent
  // acquisitions will be similar).
  ISMRMRD::Acquisition acq;
  d->readAcquisition(0, acq);
  uword nro = acq.number_of_samples();
  uword nc = acq.active_channels();

	// Handle Nx, Ny, Nz
	if(!vm.count("Nx")) {
		Nx = hdr.encoding[0].encodedSpace.matrixSize.x;
	} 
	if(!vm.count("Ny")) {
		Ny = hdr.encoding[0].encodedSpace.matrixSize.y;
	} 
	if(!vm.count("Nz")) {
		Nz = hdr.encoding[0].encodedSpace.matrixSize.z;
	} 

  Col<float> ix, iy, iz;
  Col<float> ImgCoord;
  ImgCoord = getISMRMRDImgCoord<float>(d);

  // Check and abort if we have more than one encoding space (No Navigators for
  // now).
  std::cout << "hdr.encoding.size() = " << hdr.encoding.size() << std::endl;
  if (hdr.encoding.size() != 1) {
    std::cout << "There are " << hdr.encoding.size()
              << " encoding spaces in this file" << std::endl;
    std::cout
        << "This recon does not handle more than one encoding space. Aborting."
        << std::endl;
    return -1;
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

  std::string baseFilename = "img";
  std::string filename;
  if (!outputImageFilePath.empty() && *outputImageFilePath.rbegin() != '/') {
    	outputImageFilePath += '/';
	}

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

  // Image vector for conversion to Numpy Array
  cmplx_vec img_data;
  int idx;
  std::string idx_str;

  for (uword ii = 0; ii < (*taskList)[world.rank()].size(); ii++) {
    
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
    Col<float> kx(nro), ky(nro), kz(nro), tvec(nro), k2nd_1(nro), k2nd_2(nro), k2nd_3(nro), k2nd_4(nro), k2nd_5(nro);
    Col<std::complex<float>> data(nro * nc);
    Col<std::complex<float>> ImageTemp(Nx * Ny * Nz);

                      filename = outputImageFilePath + baseFilename + "_" + "Slice" + std::to_string(NSlice) +
          								"_" + "Rep" + std::to_string(NRep) + "_" + "Avg" + std::to_string(NAvg) +
          								"_" + "Echo" + std::to_string(NEcho) + "_" + "Phase" + std::to_string(NPhase);
 
	                      senSlice = getISMRMRDCompleteSENSEMap<std::complex<float>>(d, sen, NSlice, Nx*Ny*Nz);
						            fmSlice = getISMRMRDCompleteFieldMap<float>(d, FM, NSlice, (uword) (Nx*Ny*Nz));
	                      getCompleteISMRMRDAcqData2ndorder<float>(d, acqTrack, NSlice, NRep, NAvg, NEcho, NPhase, data, kx, ky,
			                    kz, tvec, k2nd_1, k2nd_2, k2nd_3, k2nd_4, k2nd_5);
                        getISMRMRDCompleteImgCoord(ix, iy, iz, d, ImgCoord, NSlice, Nx*Ny*Nz);

	                    std::cout << "Number of elements in kx = " << kx.n_rows << std::endl;
	                    std::cout << "Number of elements in ky = " << ky.n_rows << std::endl;
	                    std::cout << "Number of elements in kz = " << kz.n_rows << std::endl;
	                    std::cout << "Number of rows in data = " << data.n_rows << std::endl;
	                    std::cout << "Number of columns in data = " << data.n_cols << std::endl;

	                    QuadPenalty<float> R(Nx, Ny, Nz, beta, dims2penalize);

                      Gdft_2ndorder<float> A(kx.n_rows, Nx*Ny*Nz,kx,ky,kz,k2nd_1,k2nd_2,k2nd_3,k2nd_4,k2nd_5,ix,iy,iz,fmSlice,tvec);
                      SENSE<float, Gdft_2ndorder<float>> Sg(A, senSlice, kx.n_rows, Nx*Ny*Nz, nc);
                      ImageTemp = reconSolve<float, SENSE<float, Gdft_2ndorder<float>>,
                            QuadPenalty<float>>(data, Sg, R, kx, ky, kz, Nx,
                            Ny, Nz, tvec, NIter);

                    if (writeNifti)
                      writeNiftiMagPhsImage<float>(filename,ImageTemp,Nx,Ny,Nz);

                    // save data for Python conversion
                    idx = NSlice * NPhaseMax * NEchoMax * NAvgMax * NRepMax + NPhase * NEchoMax * NAvgMax * NRepMax + NEcho * NAvgMax * NRepMax + NAvg * NRepMax + NRep;
                    for(int ii = 0; ii < Nx * Ny * Nz; ii++)
                        img_data.push_back(static_cast<std::complex<float>>(ImageTemp(ii)));
                    
                    // Write single numpy files
                    idx_str = std::to_string(idx);
                    cnpy::npy_save(outputImageFilePath + "images_pg_"+idx_str+".npy",&img_data[0],{Nz,Ny,Nx},"w");
                    img_data.clear();
    }

  world.barrier();

  // write single files into one numpy file
  if (world.rank() == 0){
    int Nimg = NSliceMax * NPhaseMax * NEchoMax * NAvgMax * NRepMax;
    std::string filename;
    cmplx_vec img_data_write;
    for(int k = 0; k < Nimg; k++){
      // open image
      idx_str = std::to_string(k);
      filename = outputImageFilePath + "images_pg_"+idx_str+".npy";
      cnpy::NpyArray img = cnpy::npy_load(filename);
      std::complex<float>* loaded_data = img.data<std::complex<float>>();

      // Write data to vector
      for(int j = 0; j < Nx * Ny * Nz; j++)
        img_data_write.push_back(loaded_data[j]);

      // Remove single image
      std::remove(filename.c_str());
    }
    // Write images to one numpy file
    cnpy::npy_save(outputImageFilePath + "images_pg.npy",&img_data_write[0],{NSliceMax,NPhaseMax,NEchoMax,NAvgMax,NRepMax,Nz,Ny,Nx},"w");
  }

  // Close ISMRMRD::Dataset, hdr, and acqTrack
	closeISMRMRDData(d,hdr,acqTrack);

  return 0;
}
