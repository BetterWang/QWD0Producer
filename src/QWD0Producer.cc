#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"

#include "QWD0Fitter.h"

class dso_hidden QWD0Producer final : public edm::stream::EDProducer<> {
public:
	explicit QWD0Producer(const edm::ParameterSet&);

private:
	void produce(edm::Event&, const edm::EventSetup&) override;

	QWD0Fitter theVees;
};


// Constructor
QWD0Producer::QWD0Producer(const edm::ParameterSet& iConfig) :
	theVees(iConfig, consumesCollector())
{
	produces< reco::VertexCompositeCandidateCollection >();
}


//
// Methods
//

// Producer Method
void QWD0Producer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
	using namespace edm;

	// Create auto_ptr for each collection to be stored in the Event
	std::auto_ptr< reco::VertexCompositeCandidateCollection >
		D0Candidates( new reco::VertexCompositeCandidateCollection );

	// invoke the fitter which reconstructs the vertices and fills,
	//  collections of Kshorts, Lambda0s
	theVees.fitAll(iEvent, iSetup, *D0Candidates);


	// Write the collections to the Event
	std::cout << "put D0 " << D0Candidates->size() << std::endl;
	D0Candidates->shrink_to_fit();
	iEvent.put( D0Candidates );
}

//define this as a plug-in
#include "FWCore/PluginManager/interface/ModuleDef.h"

DEFINE_FWK_MODULE(QWD0Producer);
