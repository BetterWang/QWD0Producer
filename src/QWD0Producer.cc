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
	produces< reco::VertexCompositeCandidateCollection >("KPlusPiMinus");
	produces< reco::VertexCompositeCandidateCollection >("KMinusPiPlus");
	produces< reco::VertexCompositeCandidateCollection >("KPlusPiPlus");
	produces< reco::VertexCompositeCandidateCollection >("KMinusPiMinus");
}


//
// Methods
//

// Producer Method
void QWD0Producer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
	using namespace edm;

	// Create auto_ptr for each collection to be stored in the Event K-Pi charge
	//
	std::auto_ptr< reco::VertexCompositeCandidateCollection > D0pp( new reco::VertexCompositeCandidateCollection );
	std::auto_ptr< reco::VertexCompositeCandidateCollection > D0mm( new reco::VertexCompositeCandidateCollection );
	std::auto_ptr< reco::VertexCompositeCandidateCollection > D0pm( new reco::VertexCompositeCandidateCollection );
	std::auto_ptr< reco::VertexCompositeCandidateCollection > D0mp( new reco::VertexCompositeCandidateCollection );

	// invoke the fitter which reconstructs the vertices and fills,
	//  collections of Kshorts, Lambda0s
	theVees.fitAll(iEvent, iSetup, *D0pp, *D0mm, *D0pm, *D0mp);


	// Write the collections to the Event
	D0pp->shrink_to_fit();
	D0pm->shrink_to_fit();
	D0mp->shrink_to_fit();
	D0mm->shrink_to_fit();
	std::cout << "put D0pp " << D0pp->size() << std::endl;
	std::cout << "put D0pm " << D0pm->size() << std::endl;
	std::cout << "put D0mp " << D0mp->size() << std::endl;
	std::cout << "put D0mm " << D0mm->size() << std::endl;
	iEvent.put( D0pp, "KPlusPiPlus" );
	iEvent.put( D0pm, "KPlusPiMinus" );
	iEvent.put( D0mp, "KMinusPiPlus" );
	iEvent.put( D0mm, "KMinusPiMinus" );
}

//define this as a plug-in
#include "FWCore/PluginManager/interface/ModuleDef.h"

DEFINE_FWK_MODULE(QWD0Producer);
