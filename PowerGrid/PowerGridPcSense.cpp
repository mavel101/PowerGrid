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
#include "ismrmrd/dataset.h"
#include "ismrmrd/ismrmrd.h"
#include "ismrmrd/version.h"
#include "ismrmrd/xml.h"
#include "processIsmrmrd.hpp"
#include "processNIFTI.hpp"
#include <boost/program_options.hpp>

namespace po = boost::program_options;
//using namespace PowerGrid;

int main(int argc, char **argv)
{
	std::string rawDataFilePath, outputImageFilePath;
	uword Nx, Ny, Nz, NShots = 1, type, L = 0, NIter = 10;
	double beta = 0.0;
	uword dims2penalize = 3;
	po::options_description desc("Allowed options");
	desc.add_options()("help,h", "produce help message")(
			"inputData,i", po::value<std::string>(&rawDataFilePath)->required(),
			"input ISMRMRD Raw Data file")
			/*
			("inputDataNav,-N", po::value<std::string>(&rawDataNavFilePath), "input
			ISMRMRD Navigator Raw Data")
			("outputImage,o",
			po::value<std::string>(&outputImageFilePath)->required(), "output ISMRMRD
			Image file")
			("SENSEMap,S", po::value<std::string>(&senseMapFilePath),
			 "Enable SENSE recon with the specified SENSE map in ISMRMRD image
			format")
			("FieldMap,F", po::value<std::string>(&fieldMapFilePath),
			"Enable field corrected reconstruction with the specified field map in
			ISMRMRD format")
			("Precision,P", po::value<std::string>(&precisionString),
			 "Numerical precision to use, float or double currently supported")
			 */
			("Nx,x", po::value<uword>(&Nx)->required(), "Image size in X (Required)")
			("Ny,y", po::value<uword>(&Ny)->required(),
					"Image size in Y (Required)")
			("Nz,z", po::value<uword>(&Nz)->required(), "Image size in Z (Required)")
			("NShots,s", po::value<uword>(&NShots), "Number of shots per image")
			("Beta,B", po::value<double>(&beta), "Spatial regularization penalty weight")
			("Dims2Penalize,D", po::value<uword>(&dims2penalize), "Dimensions to apply regularization to (2 or 3)")
			("CGIterations,n", po::value<uword>(&NIter),
					"Number of preconditioned conjugate gradient interations for main solver");

	po::variables_map vm;

	try {

		po::store(po::parse_command_line(argc, argv, desc), vm);
		po::notify(vm);

		if (vm.count("help")) {
			std::cout << desc << std::endl;
			return 1;
		}
		/*
		if(precisionString.compare("double") ==0) {
				typedef double PGPrecision;
		} else if(precisionString.compare("float") == 0) {
				typedef float PGPrecision;
		} else {
				typedef double PGPrecision;
				std::cout << "Did not recognize precision option. Defaulting to
		double precision." << std::endl;
		}
		*/


	}
	catch (boost::program_options::error& e) {
		std::cerr << "Error: " << e.what() << std::endl;
		std::cout << desc << std::endl;
		return 1;
	}

	arma::Col<float> FM;
	arma::Col<std::complex<float>> sen;
	arma::Col<float> PMap;
	arma::Col<float> FMap;
	arma::Col<std::complex<float>> SMap;
	ISMRMRD::Dataset* d;
	ISMRMRD::IsmrmrdHeader hdr;
	acqTracking* acqTrack;
	processISMRMRDInput<float>(rawDataFilePath, d, hdr, FM, sen, acqTrack);

	//savemat("testFM.mat", "FM", FM);
	//savemat("testSen.mat", "sen", sen);

	std::cout << "Number of elements in SENSE Map = " << sen.n_rows << std::endl;
	std::cout << "Number of elements in Field Map = " << FM.n_rows << std::endl;

	uword numAcq = d->getNumberOfAcquisitions();

	// Grab first acquisition to get parameters (We assume all subsequent
	// acquisitions will be similar).
	ISMRMRD::Acquisition acq;
	d->readAcquisition(0, acq);
	uword nro = acq.number_of_samples();
	uword nc = acq.active_channels();

	Col<float> ix, iy, iz;
	initImageSpaceCoords(ix, iy, iz, Nx, Ny, Nz);
	Col<float> kx(nro), ky(nro), kz(nro), tvec(nro);
	Col<std::complex<float>> data(nro*nc);
	Col<std::complex<float>> ImageTemp(Nx*Ny*Nz);

	// Check and abort if we have more than one encoding space (No Navigators for
	// now).
	if (hdr.encoding.size()!=1) {
		std::cout << "There are " << hdr.encoding.size()
		          << " encoding spaces in this file" << std::endl;
		std::cout
				<< "This recon does not handle more than one encoding space. Aborting."
				<< std::endl;
		return -1;
	}

	int NShotMax = hdr.encoding[0].encodingLimits.kspace_encoding_step_1->maximum;
	int NParMax = hdr.encoding[0].encodingLimits.kspace_encoding_step_2->maximum;
	int NSliceMax = hdr.encoding[0].encodingLimits.slice->maximum;
	int NSetMax = hdr.encoding[0].encodingLimits.set->maximum;
	int NRepMax = hdr.encoding[0].encodingLimits.repetition->maximum;
	int NAvgMax = hdr.encoding[0].encodingLimits.average->maximum;
	int NSegMax = hdr.encoding[0].encodingLimits.segment->maximum;
	int NEchoMax = hdr.encoding[0].encodingLimits.contrast->maximum;
	int NPhaseMax = hdr.encoding[0].encodingLimits.phase->maximum;

	std::cout << "NParMax = " << NParMax << std::endl;
	std::cout << "NShotMax = " << NShotMax << std::endl;
	std::cout << "NSliceMax = " << NSliceMax << std::endl;
	std::cout << "NSetMax = " << NSetMax << std::endl;
	std::cout << "NRepMax = " << NRepMax << std::endl;
	std::cout << "NAvgMax = " << NAvgMax << std::endl;
	std::cout << "NEchoMax = " << NEchoMax << std::endl;
	std::cout << "NPhaseMax = " << NPhaseMax << std::endl;
	std::cout << "NSegMax = " << NSegMax << std::endl;
	std::cout << "About to loop through the counters and scan the file"
	          << std::endl;
  std::string baseFilename = "pcSENSE";
	std::string filename;
	uword NSet = 0; //Set is only used for arrayed ADCs
	uword NSeg = 0;
	for (uword NPhase = 0; NPhase<=NPhaseMax; NPhase++) {
		for (uword NEcho = 0; NEcho<=NEchoMax; NEcho++) {
			for (uword NAvg = 0; NAvg<=NAvgMax; NAvg++) {
				for (uword NRep = 0; NRep<NRepMax+1; NRep++) {
					for (uword NSlice = 0; NSlice<=NSliceMax; NSlice++) {

						filename = baseFilename + "_" + "Slice" + std::to_string(NSlice) +
								"_" + "Rep" + std::to_string(NRep) + "_" + "Avg" + std::to_string(NAvg) +
								"_" + "Echo" + std::to_string(NEcho) + "_" + "Phase" + std::to_string(NPhase);

						getCompleteISMRMRDAcqData<float>(d, acqTrack, NSlice, NRep, NAvg, NEcho, NPhase, data, kx, ky,
								kz, tvec);

						PMap = getISMRMRDCompletePhaseMap<float>(d, NSlice, NSet, NRep, NAvg, NPhase, NEcho, NSeg,
								(uword) (Nx*Ny*Nz));
						SMap = getISMRMRDCompleteSENSEMap<std::complex<float>>(d, NSlice, (uword) (Nx*Ny*Nz));
						FMap = getISMRMRDCompleteFieldMap<float>(d, NSlice, (uword) (Nx*Ny*Nz));

						std::cout << "Number of elements in SMap = " << SMap.n_rows << std::endl;
						std::cout << "Number of elements in kx = " << kx.n_rows << std::endl;
						std::cout << "Number of elements in ky = " << ky.n_rows << std::endl;
						std::cout << "Number of elements in kz = " << kz.n_rows << std::endl;
						std::cout << "Number of rows in phase map = " << PMap.n_rows << std::endl;
						std::cout << "Number of rows in data = " << data.n_rows << std::endl;

						std::cout << "Number of columns in data = " << data.n_cols << std::endl;

						//Gnufft<float> A(kx.n_rows, (float) 2.0, Nx, Ny, Nz, kx, ky, kz, ix,
						// iy, iz);
						//Gdft<float> A(kx.n_rows, Nx*Ny*Nz,kx,ky,kz,ix,iy,iz,FM,tvec);
						pcSENSE<float> S_DWI(kx, ky, kz, Nx, Ny, Nz, nc, tvec, SMap, FMap,
								0-PMap);
						//pcSENSE<float, Gnufft<float>> Sg(A, sen, kx.n_rows, Nx*Ny*Nz, nc);
						QuadPenalty<float> R(Nx, Ny, Nz, beta, dims2penalize);

						ImageTemp = reconSolve<float, pcSENSE<float>, QuadPenalty<float>>(data, S_DWI, R, kx, ky, kz,
								Nx,
								Ny, Nz, tvec, NIter);
						//writeISMRMRDImageData<float>(d, ImageTemp, Nx, Ny, Nz);
						writeNiftiMagPhsImage<float>(filename,ImageTemp,Nx,Ny,Nz);
					}
				}
			}
		}
	}

  // Close ISMRMRD::Dataset
  delete d;

  return 0;
}
