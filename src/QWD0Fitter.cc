#include "QWD0Fitter.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/PatternTools/interface/TSCBLBuilderNoMaterial.h"
#include <Math/Functions.h>
#include <Math/SVector.h>
#include <Math/SMatrix.h>
#include <typeinfo>
#include <memory>
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "CommonTools/CandUtils/interface/AddFourMomenta.h"

// pdg mass constants
namespace {
	const double piMass = 0.13957018;
	const double piMassSquared = piMass*piMass;
	const double protonMass = 0.938272046;
	const double protonMassSquared = protonMass*protonMass;
	const double kaonMass = 0.493667;
	const double kaonMassSquared = kaonMass*kaonMass;
	const double kShortMass = 0.497614;
	const double lambdaMass = 1.115683;
	const double D0Mass = 1.86484;
}

typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;
typedef ROOT::Math::SVector<double, 3> SVector3;

QWD0Fitter::QWD0Fitter(const edm::ParameterSet& theParameters, edm::ConsumesCollector && iC)
{
	token_beamSpot = iC.consumes<reco::BeamSpot>(theParameters.getParameter<edm::InputTag>("beamSpot"));
	useVertex_ = theParameters.getParameter<bool>("useVertex");
	token_vertices = iC.consumes<std::vector<reco::Vertex>>(theParameters.getParameter<edm::InputTag>("vertices"));

	token_tracks = iC.consumes<reco::TrackCollection>(theParameters.getParameter<edm::InputTag>("trackRecoAlgorithm"));
	vertexFitter_ = theParameters.getParameter<bool>("vertexFitter");
	useRefTracks_ = theParameters.getParameter<bool>("useRefTracks");

	// cuts on initial track selection
	tkChi2Cut_ = theParameters.getParameter<double>("tkChi2Cut");
	tkNHitsCut_ = theParameters.getParameter<int>("tkNHitsCut");
	tkPtCut_ = theParameters.getParameter<double>("tkPtCut");
	tkIPSigXYCut_ = theParameters.getParameter<double>("tkIPSigXYCut");
	tkIPSigZCut_ = theParameters.getParameter<double>("tkIPSigZCut");

	// cuts on vertex
	vtxChi2Cut_ = theParameters.getParameter<double>("vtxChi2Cut");
	vtxDecaySigXYZCut_ = theParameters.getParameter<double>("vtxDecaySigXYZCut");
	vtxDecaySigXYCut_ = theParameters.getParameter<double>("vtxDecaySigXYCut");
	// miscellaneous cuts
	tkDCACut_ = theParameters.getParameter<double>("tkDCACut");
	mPiPiCut_ = theParameters.getParameter<double>("mPiPiCut");
	innerHitPosCut_ = theParameters.getParameter<double>("innerHitPosCut");
	cosThetaXYCut_ = theParameters.getParameter<double>("cosThetaXYCut");
	cosThetaXYZCut_ = theParameters.getParameter<double>("cosThetaXYZCut");
	// cuts on the D0 candidate mass
	D0MassCut_ = theParameters.getParameter<double>("D0MassCut");
}

// method containing the algorithm for vertex reconstruction
// K-Pi charge
void QWD0Fitter::fitAll(const edm::Event& iEvent, const edm::EventSetup& iSetup,
		reco::VertexCompositeCandidateCollection & d0s)
{
	using std::vector;

	edm::Handle<reco::TrackCollection> theTrackHandle;
	iEvent.getByToken(token_tracks, theTrackHandle);
	if (!theTrackHandle->size()) return;
	const reco::TrackCollection* theTrackCollection = theTrackHandle.product();

	edm::Handle<reco::BeamSpot> theBeamSpotHandle;
	iEvent.getByToken(token_beamSpot, theBeamSpotHandle);
	const reco::BeamSpot* theBeamSpot = theBeamSpotHandle.product();
	math::XYZPoint referencePos(theBeamSpot->position());

	reco::Vertex referenceVtx;
	if (useVertex_) {
		edm::Handle<std::vector<reco::Vertex>> vertices;
		iEvent.getByToken(token_vertices, vertices);
		referenceVtx = vertices->at(0);
		referencePos = referenceVtx.position();
	}

	edm::ESHandle<MagneticField> theMagneticFieldHandle;
	iSetup.get<IdealMagneticFieldRecord>().get(theMagneticFieldHandle);
	const MagneticField* theMagneticField = theMagneticFieldHandle.product();

	std::vector<reco::TrackRef> theTrackRefs;
	std::vector<reco::TransientTrack> theTransTracks;

	// fill vectors of TransientTracks and TrackRefs after applying preselection cuts
	for (reco::TrackCollection::const_iterator iTk = theTrackCollection->begin(); iTk != theTrackCollection->end(); ++iTk) {
		const reco::Track* tmpTrack = &(*iTk);
		double ipsigXY = std::abs(tmpTrack->dxy(*theBeamSpot)/tmpTrack->dxyError());
		if (useVertex_) ipsigXY = std::abs(tmpTrack->dxy(referencePos)/tmpTrack->dxyError());
		double ipsigZ = std::abs(tmpTrack->dz(referencePos)/tmpTrack->dzError());
		if (tmpTrack->normalizedChi2() < tkChi2Cut_ && tmpTrack->numberOfValidHits() >= tkNHitsCut_ &&
				tmpTrack->pt() > tkPtCut_ && ipsigXY > tkIPSigXYCut_ && ipsigZ > tkIPSigZCut_) {
			reco::TrackRef tmpRef(theTrackHandle, std::distance(theTrackCollection->begin(), iTk));
			theTrackRefs.push_back(std::move(tmpRef));
			reco::TransientTrack tmpTransient(*tmpRef, theMagneticField);
			theTransTracks.push_back(std::move(tmpTransient));
		}
	}
	// good tracks have now been selected for vertexing

	// loop over tracks and vertex good charged track pairs
	for (unsigned int trdx1 = 0; trdx1 < theTrackRefs.size(); ++trdx1) {
		for (unsigned int trdx2 = trdx1 + 1; trdx2 < theTrackRefs.size(); ++trdx2) {

			reco::TrackRef TrackRef1 = theTrackRefs[trdx1];
			reco::TrackRef TrackRef2 = theTrackRefs[trdx2];
			reco::TransientTrack* TransTkPtr1 = &theTransTracks[trdx1];
			reco::TransientTrack* TransTkPtr2 = &theTransTracks[trdx2];

			if (theTrackRefs[trdx1]->charge() == 0 or  theTrackRefs[trdx2]->charge() == 0) continue;
			//int chIdx = (theTrackRefs[trdx1]->charge()>0?1:0)*2 + (theTrackRefs[trdx2]->charge()>0?1:0);
			int charge1 = theTrackRefs[trdx1]->charge() > 0 ? 1:-1;
			int charge2 = theTrackRefs[trdx2]->charge() > 0 ? 1:-1;

			// measure distance between tracks at their closest approach
			if (!TransTkPtr1->impactPointTSCP().isValid() || !TransTkPtr2->impactPointTSCP().isValid()) {
				std::cout << " --> " << __LINE__ << std::endl;
				continue;
			}
			FreeTrajectoryState const & State1 = TransTkPtr1->impactPointTSCP().theState();
			FreeTrajectoryState const & State2 = TransTkPtr2->impactPointTSCP().theState();
			ClosestApproachInRPhi cApp;
			cApp.calculate(State1, State2);
			if (!cApp.status()) {
				std::cout << " --> " << __LINE__ << " cApp.status() = " << cApp.status() << std::endl;
				continue;
			}
			float dca = std::abs(cApp.distance());
			if (dca > tkDCACut_) {
				// XXX
				std::cout << " --> " << __LINE__ << " dca = " << dca << std::endl;
				continue;
			}
			// the POCA should at least be in the sensitive volume
			GlobalPoint cxPt = cApp.crossingPoint();
			if (sqrt(cxPt.x()*cxPt.x() + cxPt.y()*cxPt.y()) > 120. || std::abs(cxPt.z()) > 300.) {
				std::cout << " --> " << __LINE__ << std::endl;
				continue;
			}

			// the tracks should at least point in the same quadrant
			TrajectoryStateClosestToPoint const & TSCP1 = TransTkPtr1->trajectoryStateClosestToPoint(cxPt);
			TrajectoryStateClosestToPoint const & TSCP2 = TransTkPtr2->trajectoryStateClosestToPoint(cxPt);
			if (!TSCP1.isValid() || !TSCP2.isValid()) {
				std::cout << " --> " << __LINE__ << std::endl;
				continue;
			}
			if (TSCP1.momentum().dot(TSCP2.momentum())  < 0) {
				std::cout << " --> " << __LINE__ << std::endl;
				continue;
			}

			// calculate mPiPi
			double totalE = sqrt(TSCP1.momentum().mag2() + piMassSquared) + sqrt(TSCP2.momentum().mag2() + piMassSquared);
			double totalESq = totalE*totalE;
			double totalPSq = (TSCP1.momentum() + TSCP2.momentum()).mag2();
			double mass = sqrt(totalESq - totalPSq);
			if (mass > mPiPiCut_) {
				std::cout << " --> " << __LINE__ << " mPiPi = " << mass << std::endl;
				continue;
			}

			// Fill the vector of TransientTracks to send to KVF
			std::vector<reco::TransientTrack> transTracks;
			transTracks.reserve(2);
			transTracks.push_back(*TransTkPtr1);
			transTracks.push_back(*TransTkPtr2);

			// create the vertex fitter object and vertex the tracks
			TransientVertex theRecoVertex;
			if (vertexFitter_) {
				KalmanVertexFitter theKalmanFitter(useRefTracks_ == 0 ? false : true);
				theRecoVertex = theKalmanFitter.vertex(transTracks);
			} else if (!vertexFitter_) {
				useRefTracks_ = false;
				AdaptiveVertexFitter theAdaptiveFitter;
				theRecoVertex = theAdaptiveFitter.vertex(transTracks);
			}
			if (!theRecoVertex.isValid()) {
				std::cout << " --> " << __LINE__ << " failed KalmanVertexFitter" << std::endl;
				continue;
			}

			reco::Vertex theVtx = theRecoVertex;
			if (theVtx.normalizedChi2() > vtxChi2Cut_) {
				// XXX
				std::cout << " --> " << __LINE__ << " theVtx.normalizedChi2() " << theVtx.normalizedChi2() << std::endl;
				continue;
			}
			GlobalPoint vtxPos(theVtx.x(), theVtx.y(), theVtx.z());

			// 2D decay significance
			SMatrixSym3D totalCov = theBeamSpot->rotatedCovariance3D() + theVtx.covariance();
			if (useVertex_) totalCov = referenceVtx.covariance() + theVtx.covariance();
			SVector3 distVecXY(vtxPos.x()-referencePos.x(), vtxPos.y()-referencePos.y(), 0.);
			double distMagXY = ROOT::Math::Mag(distVecXY);
			double sigmaDistMagXY = sqrt(ROOT::Math::Similarity(totalCov, distVecXY)) / distMagXY;
			if (distMagXY/sigmaDistMagXY < vtxDecaySigXYCut_) {
				// XXX
				std::cout << " --> " << __LINE__ << " distMagXY/sigmaDistMagXY = " << distMagXY/sigmaDistMagXY << std::endl;
				continue;
			}

			// 3D decay significance
			SVector3 distVecXYZ(vtxPos.x()-referencePos.x(), vtxPos.y()-referencePos.y(), vtxPos.z()-referencePos.z());
			double distMagXYZ = ROOT::Math::Mag(distVecXYZ);
			double sigmaDistMagXYZ = sqrt(ROOT::Math::Similarity(totalCov, distVecXYZ)) / distMagXYZ;
			if (distMagXYZ/sigmaDistMagXYZ < vtxDecaySigXYZCut_) {
				// XXX
				std::cout << " --> " << __LINE__ << " distMagXYZ/sigmaDistMagXYZ = " << distMagXYZ/sigmaDistMagXYZ << std::endl;
				continue;
			}

		// make sure the vertex radius is within the inner track hit radius
		// comment out for missing TrackExtra
//			if (innerHitPosCut_ > 0. && TrackRef1->innerOk()) {
//				reco::Vertex::Point posTkHitPos = TrackRef1->innerPosition();
//				double posTkHitPosD2 =  (posTkHitPos.x()-referencePos.x())*(posTkHitPos.x()-referencePos.x()) +
//					(posTkHitPos.y()-referencePos.y())*(posTkHitPos.y()-referencePos.y());
//				if (sqrt(posTkHitPosD2) < (distMagXY - sigmaDistMagXY*innerHitPosCut_)) continue;
//			}
//			if (innerHitPosCut_ > 0. && TrackRef2->innerOk()) {
//				reco::Vertex::Point negTkHitPos = TrackRef2->innerPosition();
//				double negTkHitPosD2 = (negTkHitPos.x()-referencePos.x())*(negTkHitPos.x()-referencePos.x()) +
//					(negTkHitPos.y()-referencePos.y())*(negTkHitPos.y()-referencePos.y());
//				if (sqrt(negTkHitPosD2) < (distMagXY - sigmaDistMagXY*innerHitPosCut_)) continue;
//			}

			std::auto_ptr<TrajectoryStateClosestToPoint> traj1;
			std::auto_ptr<TrajectoryStateClosestToPoint> traj2;
			std::vector<reco::TransientTrack> theRefTracks;
			if (theRecoVertex.hasRefittedTracks()) {
				theRefTracks = theRecoVertex.refittedTracks();
			}

			if (useRefTracks_ && theRefTracks.size() > 1) {
				reco::TransientTrack* thePositiveRefTrack = &(theRefTracks[0]);
				reco::TransientTrack* theNegativeRefTrack = &(theRefTracks[1]);
				if (thePositiveRefTrack == 0 || theNegativeRefTrack == 0) continue;
				traj1.reset(new TrajectoryStateClosestToPoint(thePositiveRefTrack->trajectoryStateClosestToPoint(vtxPos)));
				traj2.reset(new TrajectoryStateClosestToPoint(theNegativeRefTrack->trajectoryStateClosestToPoint(vtxPos)));
			} else {
				traj1.reset(new TrajectoryStateClosestToPoint(TransTkPtr1->trajectoryStateClosestToPoint(vtxPos)));
				traj2.reset(new TrajectoryStateClosestToPoint(TransTkPtr2->trajectoryStateClosestToPoint(vtxPos)));
			}

			if (traj1.get() == 0 || traj2.get() == 0 || !traj1->isValid() || !traj2->isValid()) continue;

			GlobalVector P1(traj1->momentum());
			GlobalVector P2(traj2->momentum());
			GlobalVector totalP(P1 + P2);

			// 2D pointing angle
			double dx = theVtx.x()-referencePos.x();
			double dy = theVtx.y()-referencePos.y();
			double px = totalP.x();
			double py = totalP.y();
			double angleXY = (dx*px+dy*py)/(sqrt(dx*dx+dy*dy)*sqrt(px*px+py*py));
			if (angleXY < cosThetaXYCut_) {
				// XXX
				std::cout << " --> " << __LINE__ << " angleXY = " << angleXY << std::endl;
				continue;
			}

			// 3D pointing angle
			double dz = theVtx.z()-referencePos.z();
			double pz = totalP.z();
			double angleXYZ = (dx*px+dy*py+dz*pz)/(sqrt(dx*dx+dy*dy+dz*dz)*sqrt(px*px+py*py+pz*pz));
			if (angleXYZ < cosThetaXYZCut_) {
				// XXX
				std::cout << " --> " << __LINE__ << " angleXYZ = " << angleXYZ << std::endl;
				continue;
			}

			// calculate total energy of D0
			double piE1 = sqrt(P1.mag2() + piMassSquared);
			double piE2 = sqrt(P2.mag2() + piMassSquared);
			double kaonE1 = sqrt(P1.mag2() + kaonMassSquared);
			double kaonE2 = sqrt(P2.mag2() + kaonMassSquared);

			double	D0Epk = piE1 + kaonE2;
			double	D0Ekp = kaonE1 + piE2;

			// Create momentum 4-vectors for the 3 candidate types
			const reco::Particle::LorentzVector D0P4pk(totalP.x(), totalP.y(), totalP.z(), D0Epk);
			const reco::Particle::LorentzVector D0P4kp(totalP.x(), totalP.y(), totalP.z(), D0Ekp);

			reco::Particle::Point vtx(theVtx.x(), theVtx.y(), theVtx.z());
			const reco::Vertex::CovarianceMatrix vtxCov(theVtx.covariance());
			double vtxChi2(theVtx.chi2());
			double vtxNdof(theVtx.ndof());

			// Create the VertexCompositeCandidate object that will be stored in the Event
			auto theD0pk = new reco::VertexCompositeCandidate(0, D0P4pk, vtx, vtxCov, vtxChi2, vtxNdof);
			auto theD0kp = new reco::VertexCompositeCandidate(0, D0P4kp, vtx, vtxCov, vtxChi2, vtxNdof);

			// Create daughter candidates for the VertexCompositeCandidates
			reco::RecoChargedCandidate thePiCand1(charge1, reco::Particle::LorentzVector(P1.x(), P1.y(), P1.z(), piE1), vtx);
			thePiCand1.setTrack(TrackRef1);

			reco::RecoChargedCandidate thePiCand2(charge2, reco::Particle::LorentzVector(P2.x(), P2.y(), P2.z(), piE2), vtx);
			thePiCand2.setTrack(TrackRef2);

			reco::RecoChargedCandidate theKaonCand1(charge1, reco::Particle::LorentzVector(P1.x(), P1.y(), P1.z(), kaonE1), vtx);
			theKaonCand1.setTrack(TrackRef1);

			reco::RecoChargedCandidate theKaonCand2(charge2, reco::Particle::LorentzVector(P2.x(), P2.y(), P2.z(), kaonE2), vtx);
			theKaonCand2.setTrack(TrackRef2);

			AddFourMomenta addp4;
			// Store the daughter Candidates in the VertexCompositeCandidates if they pass mass cuts
			theD0pk->addDaughter(thePiCand1);
			theD0pk->addDaughter(theKaonCand2);
			addp4.set(*theD0pk);
			std::cout << " --> " << __LINE__ << " theD0pk->mass() = " << theD0pk->mass() << std::endl;
			if ( theD0pk->mass() < D0Mass + D0MassCut_ and theD0pk->mass() > D0Mass - D0MassCut_ ) {
				if ( charge1 > 0 and charge2 < 0 ) {
					theD0pk->setPdgId(421);
				} else
				if ( charge1 < 0 and charge2 > 0 ) {
					theD0pk->setPdgId(-421);
				} else
				if ( charge1 > 0 and charge2 > 0 ) {
					theD0pk->setPdgId(81);
				} else
				if ( charge1 < 0 and charge2 < 0 ) {
					theD0pk->setPdgId(-81);
				}
				d0s.push_back(std::move(*theD0pk));
			}

			theD0kp->addDaughter(theKaonCand1);
			theD0kp->addDaughter(thePiCand2);
			addp4.set(*theD0kp);
			std::cout << " --> " << __LINE__ << " theD0kp->mass() = " << theD0kp->mass() << std::endl;
			if ( theD0kp->mass() < D0Mass + D0MassCut_ and theD0kp->mass() > D0Mass - D0MassCut_ ) {
				if ( charge1 > 0 and charge2 < 0 ) {
					theD0kp->setPdgId(-421);
				} else
				if ( charge1 < 0 and charge2 > 0 ) {
					theD0kp->setPdgId(421);
				} else
				if ( charge1 > 0 and charge2 > 0 ) {
					theD0kp->setPdgId(81);
				} else
				if ( charge1 < 0 and charge2 < 0 ) {
					theD0kp->setPdgId(-81);
				}
				d0s.push_back(std::move(*theD0kp));
			}
			delete theD0pk;
			delete theD0kp;
			theD0pk = theD0kp = nullptr;
		}
	}
}

