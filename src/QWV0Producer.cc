#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
//#include "DataFormats/V0Candidate/interface/V0Candidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"

#include "QWV0Fitter.h"

class dso_hidden QWV0Producer final : public edm::stream::EDProducer<> {
public:
	explicit QWV0Producer(const edm::ParameterSet&);

private:
	void produce(edm::Event&, const edm::EventSetup&) override;

	bool doKShorts_;
	bool doLambdas_;
	bool doD0s_;

	V0Fitter theVees;
};


// Constructor
QWV0Producer::QWV0Producer(const edm::ParameterSet& iConfig) :
	theVees(iConfig, consumesCollector())
{
	doKShorts_ = iConfig.getParameter<bool>("doKShorts");
	doLambdas_ = iConfig.getParameter<bool>("doLambdas");
	doD0s_ = iConfig.getParameter<bool>("doD0s");
	if ( doKShorts_ ) produces< reco::VertexCompositeCandidateCollection >("Kshort");
	if ( doLambdas_ ) produces< reco::VertexCompositeCandidateCollection >("Lambda");
	if ( doD0s_ ) produces< reco::VertexCompositeCandidateCollection >("D0");
}


//
// Methods
//

// Producer Method
void QWV0Producer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
	using namespace edm;

	// Create auto_ptr for each collection to be stored in the Event
	std::auto_ptr< reco::VertexCompositeCandidateCollection > 
		kShortCandidates( new reco::VertexCompositeCandidateCollection );
	std::auto_ptr< reco::VertexCompositeCandidateCollection >
		lambdaCandidates( new reco::VertexCompositeCandidateCollection );
	std::auto_ptr< reco::VertexCompositeCandidateCollection >
		D0Candidates( new reco::VertexCompositeCandidateCollection );

	// invoke the fitter which reconstructs the vertices and fills,
	//  collections of Kshorts, Lambda0s
	theVees.fitAll(iEvent, iSetup, *kShortCandidates, *lambdaCandidates, *D0Candidates);


	// Write the collections to the Event
	if ( doKShorts_ ) {
		kShortCandidates->shrink_to_fit();
		iEvent.put( kShortCandidates, std::string("Kshort") );
	}
	if ( doLambdas_ ) {
		lambdaCandidates->shrink_to_fit();
		iEvent.put( lambdaCandidates, std::string("Lambda") );
	}
	if ( doD0s_ ) {
		D0Candidates->shrink_to_fit();
		iEvent.put( lambdaCandidates, std::string("D0") );
	}

}

//define this as a plug-in
#include "FWCore/PluginManager/interface/ModuleDef.h"

DEFINE_FWK_MODULE(QWV0Producer);
//DEFINE_FWK_MODULE(V0finder);
