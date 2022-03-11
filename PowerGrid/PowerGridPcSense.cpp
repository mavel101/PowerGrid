/*
(C) Copyright 2015-2016 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/*****************************************************************************

    File Name   [PowerGridPcSense.cpp]

    Synopsis    [PowerGrid reconstruction executable supporting ISMRMRD format
                                as input and Phase Corrected Sense Algorithm.]

    Description [This reconstruction supports 2D and 3D reconstructions with
                                time segmentation for field correction. This
 support is experimental.]

    Revision    [0.1.0; Alex Cerjanic, BIOE UIUC]

    Date        [06/19/2017]

 *****************************************************************************/

// //Project headers.
#include "PowerGrid.h"
#include "processIsmrmrd.hpp"
#include "processNIFTI.hpp"
#include <boost/program_options.hpp>
#include <chrono>
#include "cnpy.h"

namespace po = boost::program_options;
//using namespace PowerGrid;

typedef std::vector<std::complex<float>> cmplx_vec;

int main(int argc, char **argv) {
	std::string rawDataFilePath, outputImageFilePath;
	uword Nx, Ny, Nz, NShots = 1, NIter = 10;
	//uword type, L = 0;
	double beta = 0.0;
	uword dims2penalize = 3;
	bool writeNifti;
	po::options_description desc("Allowed options");

	desc.add_options()("help,h", "produce help message")
			("inputData,i", po::value<std::string>(&rawDataFilePath)->required(),
			"input ISMRMRD Raw Data file")
 			("outputImage,o", po::value<std::string>(&outputImageFilePath)->required(), "output file path for Numpy files and NIFTIimages")
      		("writeNifti,w", po::bool_switch(&writeNifti)->default_value(false), "Write each image in own Nifti file if selected.")
			("Nx,x", po::value<uword>(&Nx), "Image size in X")
			("Ny,y", po::value<uword>(&Ny), "Image size in Y")
			("Nz,z", po::value<uword>(&Nz), "Image size in Z")
			("NShots,s", po::value<uword>(&NShots), "Number of shots per image")
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
	arma::Col<float> PMap;

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
		return -1;
	}

	int NShotMax  = hdr.encoding[0].encodingLimits.kspace_encoding_step_1->maximum;
	int NParMax   = hdr.encoding[0].encodingLimits.kspace_encoding_step_2->maximum;
	int NSliceMax = hdr.encoding[0].encodingLimits.slice->maximum;
	int NSetMax   = hdr.encoding[0].encodingLimits.set->maximum;
	int NRepMax   = hdr.encoding[0].encodingLimits.repetition->maximum;
	int NAvgMax   = hdr.encoding[0].encodingLimits.average->maximum;
	int NSegMax   = hdr.encoding[0].encodingLimits.segment->maximum;
	int NEchoMax  = hdr.encoding[0].encodingLimits.contrast->maximum;
	int NPhaseMax = hdr.encoding[0].encodingLimits.phase->maximum;

	std::cout << "NParMax = "   << NParMax << std::endl;
	std::cout << "NShotMax = "  << NShotMax << std::endl;
	std::cout << "NSliceMax = " << NSliceMax << std::endl;
	std::cout << "NSetMax = "   << NSetMax << std::endl;
	std::cout << "NRepMax = "   << NRepMax << std::endl;
	std::cout << "NAvgMax = "   << NAvgMax << std::endl;
	std::cout << "NEchoMax = "  << NEchoMax << std::endl;
	std::cout << "NPhaseMax = " << NPhaseMax << std::endl;
	std::cout << "NSegMax = "   << NSegMax << std::endl;

	std::cout << "About to loop through the counters and scan the file"
	          << std::endl;

  	std::string baseFilename = "img";


	std::string filename;
	if (!outputImageFilePath.empty() && *outputImageFilePath.rbegin() != '/') {
    	outputImageFilePath += '/';
	}
	
	// Image vector for conversion to Numpy Array
  	cmplx_vec img_data;

	uword NSet = 0; //Set is only used for arrayed ADCs
	uword NSeg = 0;

	for (uword NSlice = 0; NSlice<=NSliceMax; NSlice++) {


	Col<float> fmSlice;
 	Col<std::complex<float>> senSlice;
    Col<float> kx(nro), ky(nro), kz(nro), tvec(nro);
    Col<std::complex<float>> data(nro * nc);
    Col<std::complex<float>> ImageTemp(Nx * Ny * Nz);

	for (uword NPhase = 0; NPhase <= NPhaseMax; NPhase++) {
	  for (uword NEcho = 0; NEcho <= NEchoMax; NEcho++) {
        for (uword NAvg = 0; NAvg <= NAvgMax; NAvg++) {
          for (uword NRep = 0; NRep < NRepMax +1; NRep++) {


			filename = outputImageFilePath + baseFilename + "_" + "Slice" + std::to_string(NSlice) +
					"_" + "Rep" + std::to_string(NRep) + "_" + "Avg" + std::to_string(NAvg) +
					"_" + "Echo" + std::to_string(NEcho) + "_" + "Phase" + std::to_string(NPhase);

			senSlice = getISMRMRDCompleteSENSEMap<std::complex<float>>(d, sen, NSlice, Nx*Ny*Nz);
			fmSlice = getISMRMRDCompleteFieldMap<float>(d, FM, NSlice, (uword) (Nx*Ny*Nz));
			PMap = getISMRMRDCompletePhaseMap<float>(d, NSlice, NSet, NRep, NAvg, NPhase, NEcho, NSeg,
					(uword) (Nx*Ny*Nz));
			getCompleteISMRMRDAcqData<float>(d, acqTrack, NSlice, NRep, NAvg, NEcho, NPhase, data, kx, ky,
					kz, tvec);

			std::cout << "Number of elements in kx = " << kx.n_rows << std::endl;
			std::cout << "Number of elements in ky = " << ky.n_rows << std::endl;
			std::cout << "Number of elements in kz = " << kz.n_rows << std::endl;
			std::cout << "Number of rows in phase map = " << PMap.n_rows << std::endl;
			std::cout << "Number of rows in data = " << data.n_rows << std::endl;
			std::cout << "Number of columns in data = " << data.n_cols << std::endl;

			pcSENSE<float> S_DWI(kx, ky, kz, Nx, Ny, Nz, nc, tvec, senSlice, fmSlice,
					0-PMap);
			QuadPenalty<float> R(Nx, Ny, Nz, beta, dims2penalize);

			ImageTemp = reconSolve<float, pcSENSE<float>, QuadPenalty<float>>(data, S_DWI, R, kx, ky, kz,
					Nx,	Ny, Nz, tvec, NIter);
			
			if (writeNifti)
				writeNiftiMagPhsImage<float>(filename,ImageTemp,Nx,Ny,Nz);

			// save data for Python conversion
			for(int ii = 0; ii < Nx * Ny * Nz; ii++)
				img_data.push_back(static_cast<std::complex<float>>(ImageTemp(ii)));  

			}

	      }
		}
	  }
	}

  cnpy::npy_save(outputImageFilePath + "images_pg.npy",&img_data[0],{NSliceMax+1,NPhaseMax+1,NEchoMax+1,NAvgMax+1,NRepMax+1,Nz,Ny,Nx},"w");

  // Close ISMRMRD::Dataset, hdr, and acqTrack
  closeISMRMRDData(d,hdr,acqTrack);

  return 0;
}
