// system include files
#include <iostream>
#include <memory>
#include "TMath.h"
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"  

#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "LOTable.h"
#include "jetResolution.h"
#include "getJesUnc.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"





#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include<algorithm>
#define Pi 3.141593


struct sortPt
{
   bool operator()(TLorentzVector* s1, TLorentzVector* s2) const
   {
      return s1->Pt() >= s2->Pt();
   }
} mysortPt;

//
// class declaration
//

class PKUTreeMaker : public edm::EDAnalyzer {
public:
  explicit PKUTreeMaker(const edm::ParameterSet&);
  ~PKUTreeMaker();
  //static void fillDescriptions(edm::ConfigurationDescriptions & descriptions);
  
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
 // virtual void addTypeICorr( edm::Event const & event );
 // virtual double getJEC( reco::Candidate::LorentzVector& rawJetP4, const pat::Jet& jet, double& jetCorrEtaMax, std::vector<std::string> jecPayloadNames_ );
//  virtual double getJECOffset( reco::Candidate::LorentzVector& rawJetP4, const pat::Jet& jet, double& jetCorrEtaMax, std::vector<std::string> jecPayloadNames_ );
    
  bool hasMatchedPromptElectron(const reco::SuperClusterRef &sc, const edm::Handle<edm::View<pat::Electron> > &eleCol,
                                  const edm::Handle<reco::ConversionCollection> &convCol, const math::XYZPoint &beamspot,
                                  float lxyMin=2.0, float probMin=1e-6, unsigned int nHitsBeforeVtxMax=0);
  float EAch(float x); 
  float EAnh(float x);
  float EApho(float x);
    
  //  std::vector<std::string>                    jecAK4PayloadNames_;
   // boost::shared_ptr<FactorizedJetCorrector>   jecAK4_            ;
    std::vector<std::string> offsetCorrLabel_;
    
    FactorizedJetCorrector* jecOffset_;

    
    edm::Handle< double >  rho_;
//    edm::InputTag  METsRawLabel_;

    //edm::EDGetTokenT<pat::METCollection>  metInputToken_;
    //edm::EDGetTokenT<pat::METCollection>  reclusteredmetInputToken_;
    //std::vector<edm::EDGetTokenT<pat::METCollection>> mettokens;
    
    //edm::EDGetTokenT<pat::METCollection> reclusteredmetToken_;
    //edm::EDGetTokenT<pat::METCollection> metToken_;
    edm::EDGetTokenT<double> rhoToken_;

    edm::EDGetTokenT<edm::View<pat::Electron> > electronToken_ ;
    edm::EDGetTokenT<edm::View<pat::Photon> > photonToken_;
    edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
    edm::EDGetTokenT<std::vector<reco::Conversion> > conversionsToken_;

    
    std::vector<std::string> jetCorrLabel_;
  //  std::vector<std::string> jecAK4Labels;
   // std::vector<std::string> jecAK8Labels;

  // ----------member data ---------------------------
  TTree* outTree_;

  int nevent, run, ls;
  int nVtx;
  double triggerWeight, lumiWeight, pileupWeight;

  double ptel1, etaSC1, dEtaIn1, dPhiIn1,hOverE1;
  double full5x5_sigma1,ooEmooP1,d01,dz1,relIso1;
  int missingHits1, passConVeto1;

  double ptVlep, yVlep, phiVlep, massVlep;
  double met, metPhi, j1metPhi, j2metPhi, mtVlep;
    //Met JEC
    double METraw_et, METraw_phi, METraw_sumEt;
    double MET_et, MET_phi, MET_sumEt, MET_corrPx, MET_corrPy;
  double ptlep1, etalep1, philep1;
  int  lep;
    
    double useless;
    // AK4 Jets
    double ak4jet_pt[6],ak4jet_eta[6],ak4jet_phi[6],ak4jet_e[6];
    double ak4jet_csv[6],ak4jet_icsv[6];
    double drjetlep[6], drjetphoton[6];
    // AK4 Jets jes up
    double ak4jet_jesup_pt_[6],ak4jet_jesup_e_[6];
    // AK4 Jets jes down
    double ak4jet_jesdown_pt_[6],ak4jet_jesdown_e_[6];
    // AK4 Jets jer up
    double ak4jet_jerup_pt_[6],ak4jet_jerup_e_[6];
    // AK4 Jets jer down
    double ak4jet_jerdown_pt_[6],ak4jet_jerdown_e_[6];
    

  //Photon
  double photon_pt[6],photon_eta[6],photon_phi[6],photon_e[6];
  double drphotonlep[6];
  double photonet, photoneta, photonphi, photone;
  int iphoton;
  double drla;
  bool passEleVeto;
  double jet1pt, jet1eta, jet1phi, jet1e, jet1csv, jet1icsv;
  double jet2pt, jet2eta, jet2phi, jet2e, jet2csv, jet2icsv;
  double drj1a, drj2a, drj1l, drj2l;
  double Mjj, deltaeta, zepp;
    
    double jet1_jesup_pt, jet1_jesup_eta, jet1_jesup_phi, jet1_jesup_e;
    double jet2_jesup_pt, jet2_jesup_eta, jet2_jesup_phi, jet2_jesup_e;
    double jet_jesup_mass;
    
    double jet1_jesdown_pt, jet1_jesdown_eta, jet1_jesdown_phi, jet1_jesdown_e;
    double jet2_jesdown_pt, jet2_jesdown_eta, jet2_jesdown_phi, jet2_jesdown_e;
    double jet_jesdown_mass;
    
    double jet1_jerup_pt, jet1_jerup_eta, jet1_jerup_phi, jet1_jerup_e;
    double jet2_jerup_pt, jet2_jerup_eta, jet2_jerup_phi, jet2_jerup_e;
    double jet_jerup_mass;
    
    double jet1_jerdown_pt, jet1_jerdown_eta, jet1_jerdown_phi, jet1_jerdown_e;
    double jet2_jerdown_pt, jet2_jerdown_eta, jet2_jerdown_phi, jet2_jerdown_e;
    double jet_jerdown_mass;


  // Electron ID
  int el1passID, el1tightID,  el1mediumID, el1looseID;
  int mu1tightID;

  edm::InputTag electronIdTag_;
  void setDummyValues();

  /// Parameters to steer the treeDumper
  int originalNEvents_;
  double crossSectionPb_;
  double targetLumiInvPb_;
  std::string PKUChannel_;
  bool isGen_;
  std::string leptonicVSrc_;
  std::string ak4jetsSrc_;
  std::vector<std::string> jecAK4Labels_;
  //correction jet
  FactorizedJetCorrector* jecAK4_;
    //jes unc
  getJesUnc* jesAK4_;
    //jer unc
  LOTable* jetresolutiontable_;
  jetResolution* jerAK4_;
    
  std::string photonSrc_;
  std::string gravitonSrc_, metSrc_;
    
    std::map<std::string,double>  TypeICorrMap_;
    edm::InputTag mets_;
    edm::InputTag reclusteredmets_;
    edm::InputTag pfmets_;

};




float PKUTreeMaker::EAch( float x){
 float EA = 0.013;
 if(x>1.0)   EA = 0.0096;
 if(x>1.479) EA = 0.0107;
 if(x>2.0)   EA = 0.0077;
 if(x>2.2)   EA = 0.0088;
 if(x>2.3)   EA = 0.0065;
 if(x>2.4)   EA = 0.0030;
 return EA;
}

float PKUTreeMaker::EAnh( float x){
 float EA = 0.0056;
 if(x>1.0)   EA = 0.0107;
 if(x>1.479) EA = 0.0019;
 if(x>2.0)   EA = 0.0011;
 if(x>2.2)   EA = 0.0077;
 if(x>2.3)   EA = 0.0178;
 if(x>2.4)   EA = 0.1675;
 return EA;
}

float PKUTreeMaker::EApho( float x){
 float EA = 0.0896;
 if(x>1.0)   EA = 0.0762;
 if(x>1.479) EA = 0.0383;
 if(x>2.0)   EA = 0.0534;
 if(x>2.2)   EA = 0.0846;
 if(x>2.3)   EA = 0.1032;
 if(x>2.4)   EA = 0.1598;
 return EA;
}






//
// constructors and destructor
//
PKUTreeMaker::PKUTreeMaker(const edm::ParameterSet& iConfig)

{
  originalNEvents_ = iConfig.getParameter<int>("originalNEvents");
  crossSectionPb_  = iConfig.getParameter<double>("crossSectionPb");
  targetLumiInvPb_ = iConfig.getParameter<double>("targetLumiInvPb");
  PKUChannel_     = iConfig.getParameter<std::string>("PKUChannel");
  isGen_           = iConfig.getParameter<bool>("isGen");
  leptonicVSrc_ = iConfig.getParameter<std::string>("leptonicVSrc");
  ak4jetsSrc_      = iConfig.getParameter<std::string>("ak4jetsSrc");
  jecAK4Labels_   =  iConfig.getParameter<std::vector<std::string>>("jecAK4chsPayloadNames");
  photonSrc_      = iConfig.getParameter<std::string>("photonSrc");
  metSrc_          = iConfig.getParameter<std::string>("metSrc");
  electronIdTag_   = iConfig.getParameter<edm::InputTag>("electronIDs");
    
    
    rhoToken_  = consumes<double>(iConfig.getParameter<edm::InputTag>("rho"));
  //  metToken_ = consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("metSrc"));
    //reclusteredmetToken_ = consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("reclusteredmets"));
   // METsRawLabel_ = iConfig.getParameter<edm::InputTag>("pfmets");
    //  mettokens.push_back( metToken_ );
    //mettokens.push_back( reclusteredmetToken_ );
    
    //   metInputToken_ = mettokens[0];
    //   reclusteredmetInputToken_ = mettokens[1];
   
   //electron photon veto
   electronToken_    = (consumes<edm::View<pat::Electron> > (iConfig.getParameter<edm::InputTag>("electrons")))            ;
   photonToken_      = (consumes<edm::View<pat::Photon> >(iConfig.getParameter<edm::InputTag>("photons")))                ;
   beamSpotToken_    = (consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"))) ;
   conversionsToken_ = (consumes<std::vector<reco::Conversion> >(iConfig.getParameter<edm::InputTag>("conversions"))) ;
    

    
    
    jetCorrLabel_ = jecAK4Labels_;
    offsetCorrLabel_.push_back(jetCorrLabel_[0]);




  //now do what ever initialization is needed
  edm::Service<TFileService> fs;
  outTree_ = fs->make<TTree>("PKUCandidates","PKU Candidates");

  /// Basic event quantities
  outTree_->Branch("event"           ,&nevent         ,"event/I"          );
  outTree_->Branch("nVtx"            ,&nVtx           ,"nVtx/I"           );
  outTree_->Branch("ptVlep"          ,&ptVlep         ,"ptVlep/D"         );
  outTree_->Branch("yVlep"           ,&yVlep          ,"yVlep/D"          );
  outTree_->Branch("phiVlep"         ,&phiVlep        ,"phiVlep/D"        );
  outTree_->Branch("massVlep"        ,&massVlep       ,"massVlep/D"       );
  outTree_->Branch("mtVlep"          ,&mtVlep         ,"mtVlep/D"         );
  outTree_->Branch("lep"             ,&lep            ,"lep/I"            );

  /// AK4 Jets Info
/*  outTree_->Branch("ak4jet_pt"        , ak4jet_pt       ,"ak4jet_pt[6]/D"       );
  outTree_->Branch("ak4jet_eta"        , ak4jet_eta       ,"ak4jet_eta[6]/D"       );
  outTree_->Branch("ak4jet_phi"        , ak4jet_phi       ,"ak4jet_phi[6]/D"       );
  outTree_->Branch("ak4jet_e"        , ak4jet_e       ,"ak4jet_e[6]/D"       );
  outTree_->Branch("ak4jet_csv"        , ak4jet_csv       ,"ak4jet_csv[6]/D"       );
  outTree_->Branch("ak4jet_icsv"        , ak4jet_icsv       ,"ak4jet_icsv[6]/D"       );
    
    /// AK4 Jets jes up Info
    outTree_->Branch("ak4jet_jesup_pt"        , ak4jet_jesup_pt       ,"ak4jet_jesup_pt[6]/D"       );
    outTree_->Branch("ak4jet_jesup_e"        , ak4jet_jesup_e       ,"ak4jet_jesup_e[6]/D"       );
    
    /// AK4 Jets jes down Info
    outTree_->Branch("ak4jet_jesdown_pt"        , ak4jet_jesdown_pt       ,"ak4jet_jesdown_pt[6]/D"       );
    outTree_->Branch("ak4jet_jesdown_e"        , ak4jet_jesdown_e       ,"ak4jet_jesdown_e[6]/D"       );
    
    /// AK4 Jets jer up Info
    outTree_->Branch("ak4jet_jerup_pt"        , ak4jet_jerup_pt       ,"ak4jet_jerup_pt[6]/D"       );
    outTree_->Branch("ak4jet_jerup_e"        , ak4jet_jerup_e       ,"ak4jet_jerup_e[6]/D"       );
    
    /// AK4 Jets jer down Info
    outTree_->Branch("ak4jet_jerdown_pt"        , ak4jet_jerdown_pt       ,"ak4jet_jerdown_pt[6]/D"       );
    outTree_->Branch("ak4jet_jerdown_e"        , ak4jet_jerdown_e       ,"ak4jet_jerdown_e[6]/D"       );

*/
  /// Photon
  outTree_->Branch("photon_pt"        , photon_pt       ,"photon_pt[6]/D"       );
  outTree_->Branch("photon_eta"        , photon_eta       ,"photon_eta[6]/D"       );
  outTree_->Branch("photon_phi"        , photon_phi       ,"photon_phi[6]/D"       );
  outTree_->Branch("photon_e"        , photon_e       ,"photon_e[6]/D"       );
    outTree_->Branch("passEleVeto"        , &passEleVeto       ,"passEleVeto/O"       );
  outTree_->Branch("photonet"          ,&photonet         ,"photonet/D"         );
  outTree_->Branch("photoneta"          ,&photoneta         ,"photoneta/D"         );
  outTree_->Branch("photonphi"          ,&photonphi         ,"photonphi/D"         );
  outTree_->Branch("photone"          ,&photone         ,"photone/D"         );
  outTree_->Branch("iphoton"             ,&iphoton            ,"iphoton/I"            );
  outTree_->Branch("drla"          ,&drla         ,"drla/D"         );

  outTree_->Branch("jet1pt"          ,&jet1pt         ,"jet1pt/D"         );
  outTree_->Branch("jet1eta"          ,&jet1eta         ,"jet1eta/D"         );
  outTree_->Branch("jet1phi"          ,&jet1phi         ,"jet1phi/D"         );
  outTree_->Branch("jet1e"          ,&jet1e         ,"jet1e/D"         );
  outTree_->Branch("jet1csv"          ,&jet1csv         ,"jet1csv/D"         );
  outTree_->Branch("jet1icsv"          ,&jet1icsv         ,"jet1icsv/D"         );
  outTree_->Branch("jet2pt"          ,&jet2pt         ,"jet2pt/D"         );
  outTree_->Branch("jet2eta"          ,&jet2eta         ,"jet2eta/D"         );
  outTree_->Branch("jet2phi"          ,&jet2phi         ,"jet2phi/D"         );
  outTree_->Branch("jet2e"          ,&jet2e         ,"jet2e/D"         );
  outTree_->Branch("jet2csv"          ,&jet2csv         ,"jet2csv/D"         );
  outTree_->Branch("jet2icsv"          ,&jet2icsv         ,"jet2icsv/D"         );

  outTree_->Branch("drj1a"          ,&drj1a         ,"drj1a/D"         );
  outTree_->Branch("drj2a"          ,&drj2a         ,"drj2a/D"         );
  outTree_->Branch("drj1l"          ,&drj1l         ,"drj1l/D"         );
  outTree_->Branch("drj2l"          ,&drj2l         ,"drj2l/D"         );

  outTree_->Branch("Mjj"          ,&Mjj         ,"Mjj/D"         );
  outTree_->Branch("deltaeta"          ,&deltaeta         ,"deltaeta/D"         );
  outTree_->Branch("zepp"          ,&zepp         ,"zepp/D"         );


    outTree_->Branch("jet1_jesup_pt"          ,&jet1_jesup_pt         ,"jet1_jesup_pt/D"         );
    outTree_->Branch("jet1_jesup_e"          ,&jet1_jesup_e         ,"jet1_jesup_e/D"         );
    outTree_->Branch("jet2_jesup_pt"          ,&jet2_jesup_pt         ,"jet2_jesup_pt/D"         );
    outTree_->Branch("jet2_jesup_e"          ,&jet2_jesup_e         ,"jet2_jesup_e/D"         );
    outTree_->Branch("jet_jesup_mass"          ,&jet_jesup_mass         ,"jet_jesup_mass/D"         );
    
    outTree_->Branch("jet1_jesdown_pt"          ,&jet1_jesdown_pt         ,"jet1_jesdown_pt/D"         );
    outTree_->Branch("jet1_jesdown_e"          ,&jet1_jesdown_e         ,"jet1_jesdown_e/D"         );
    outTree_->Branch("jet2_jesdown_pt"          ,&jet2_jesdown_pt         ,"jet2_jesdown_pt/D"         );
    outTree_->Branch("jet2_jesdown_e"          ,&jet2_jesdown_e         ,"jet2_jesdown_e/D"         );
    outTree_->Branch("jet_jesdown_mass"          ,&jet_jesdown_mass         ,"jet_jesdown_mass/D"         );
    
    outTree_->Branch("jet1_jerup_pt"          ,&jet1_jerup_pt         ,"jet1_jerup_pt/D"         );
    outTree_->Branch("jet1_jerup_e"          ,&jet1_jerup_e         ,"jet1_jerup_e/D"         );
    outTree_->Branch("jet2_jerup_pt"          ,&jet2_jerup_pt         ,"jet2_jerup_pt/D"         );
    outTree_->Branch("jet2_jerup_e"          ,&jet2_jerup_e         ,"jet2_jerup_e/D"         );
    outTree_->Branch("jet_jerup_mass"          ,&jet_jerup_mass         ,"jet_jerup_mass/D"         );
    
    outTree_->Branch("jet1_jerdown_pt"          ,&jet1_jerdown_pt         ,"jet1_jerdown_pt/D"         );
    outTree_->Branch("jet1_jerdown_e"          ,&jet1_jerdown_e         ,"jet1_jerdown_e/D"         );
    outTree_->Branch("jet2_jerdown_pt"          ,&jet2_jerdown_pt         ,"jet2_jerdown_pt/D"         );
    outTree_->Branch("jet2_jerdown_e"          ,&jet2_jerdown_e         ,"jet2_jerdown_e/D"         );
    outTree_->Branch("jet_jerdown_mass"          ,&jet_jerdown_mass         ,"jet_jerdown_mass/D"         );
    

    
  /// Electron ID quantities
  outTree_->Branch("el1passID"       ,&el1passID      ,"el1passID/I"      );
  outTree_->Branch("el1tightID"       ,&el1tightID      ,"el1tightID/I"      );
  outTree_->Branch("el1mediumID"       ,&el1mediumID      ,"el1mediumID/I"      );
  outTree_->Branch("el1looseID"       ,&el1looseID      ,"el1looseID/I"      );
  outTree_->Branch("mu1tightID"       ,&mu1tightID      ,"mu1tightID/I"      );

  /// Generic kinematic quantities
  outTree_->Branch("ptlep1"          ,&ptlep1         ,"ptlep1/D"         );
  outTree_->Branch("etalep1"         ,&etalep1        ,"etalep1/D"        );
  outTree_->Branch("philep1"         ,&philep1        ,"philep1/D"        );
  outTree_->Branch("met"             ,&met            ,"met/D"            );
  outTree_->Branch("metPhi"          ,&metPhi         ,"metPhi/D"         );
  outTree_->Branch("j1metPhi"          ,&j1metPhi         ,"j1metPhi/D"         );
  outTree_->Branch("j2metPhi"          ,&j2metPhi         ,"j2metPhi/D"         );

    
    outTree_->Branch("METraw_et",&METraw_et,"METraw_et/D");
    outTree_->Branch("METraw_phi",&METraw_phi,"METraw_phi/D");
    outTree_->Branch("METraw_sumEt",&METraw_sumEt,"METraw_sumEt/D");
    outTree_->Branch("MET_et",&MET_et,"MET_et/D");
    outTree_->Branch("MET_phi",&MET_phi,"MET_phi/D");
    outTree_->Branch("MET_sumEt",&MET_sumEt,"MET_sumEt/D");
    outTree_->Branch("MET_corrPx",&MET_corrPx,"MET_corrPx/D");
    outTree_->Branch("MET_corrPy",&MET_corrPy,"MET_corrPy/D");


  /// Other quantities
  outTree_->Branch("triggerWeight"   ,&triggerWeight  ,"triggerWeight/D"  );
  outTree_->Branch("lumiWeight"      ,&lumiWeight     ,"lumiWeight/D"     );
  outTree_->Branch("pileupWeight"    ,&pileupWeight   ,"pileupWeight/D"   );
    
    
    
}


/*double PKUTreeMaker::getJEC( reco::Candidate::LorentzVector& rawJetP4, const pat::Jet& jet, double& jetCorrEtaMax, std::vector<std::string> jecPayloadNames_ ){
    
    std::vector<JetCorrectorParameters> vPar;
    //         vPar.clear();
    for ( std::vector<std::string>::const_iterator payloadBegin = jecAK4Labels_.begin(), payloadEnd = jecAK4Labels_.end(), ipayload = payloadBegin; ipayload != payloadEnd; ++ipayload ) {
        JetCorrectorParameters pars(*ipayload);
        vPar.push_back(pars);
    }
    
    jecAK4_ = new FactorizedJetCorrector(vPar);
    
    double jetCorrFactor = 1.;
    if ( fabs(rawJetP4.eta()) < jetCorrEtaMax ){
        jecAK4_->setJetEta( rawJetP4.eta() );
        jecAK4_->setJetPt ( rawJetP4.pt() );
        jecAK4_->setJetE  ( rawJetP4.energy() );
        jecAK4_->setJetPhi( rawJetP4.phi()    );
        jecAK4_->setJetA  ( jet.jetArea() );
        jecAK4_->setRho   ( *(rho_.product()) );
        jecAK4_->setNPV   ( nVtx );
        jetCorrFactor = jecAK4_->getCorrection();
    }
    
    reco::Candidate::LorentzVector corrJetP4 = rawJetP4;
    corrJetP4 *= jetCorrFactor;
    
    return jetCorrFactor;
    
}

double PKUTreeMaker::getJECOffset( reco::Candidate::LorentzVector& rawJetP4, const pat::Jet& jet, double& jetCorrEtaMax, std::vector<std::string> jecPayloadNames_ ){
    
    std::vector<JetCorrectorParameters> vPar;
    //         vPar.clear();
    for ( std::vector<std::string>::const_iterator payloadBegin = offsetCorrLabel_.begin(), payloadEnd = offsetCorrLabel_.end(), ipayload = payloadBegin; ipayload != payloadEnd; ++ipayload ) {
        JetCorrectorParameters pars(*ipayload);
        vPar.push_back(pars);
    }
    
    jecOffset_ = new FactorizedJetCorrector(vPar);
    
    double jetCorrFactor = 1.;
    if ( fabs(rawJetP4.eta()) < jetCorrEtaMax ){
        jecOffset_->setJetEta( rawJetP4.eta()     );
        jecOffset_->setJetPt ( rawJetP4.pt()      );
        jecOffset_->setJetE  ( rawJetP4.energy()  );
        jecOffset_->setJetPhi( rawJetP4.phi()     );
        jecOffset_->setJetA  ( jet.jetArea()      );
        jecOffset_->setRho   ( *(rho_.product())  );
        jecOffset_->setNPV   ( nVtx  );
        jetCorrFactor = jecOffset_->getCorrection();
    }
    
    reco::Candidate::LorentzVector corrJetP4 = rawJetP4;
    corrJetP4 *= jetCorrFactor;
    
    return jetCorrFactor;
    
}


void PKUTreeMaker::addTypeICorr( edm::Event const & event ){
    TypeICorrMap_.clear();
    edm::Handle<pat::JetCollection> jets_;
    event.getByLabel("slimmedJets", jets_);
    event.getByToken(rhoToken_      , rho_     );
    edm::Handle<reco::VertexCollection> vertices_;
    event.getByLabel("offlineSlimmedPrimaryVertices", vertices_);
    edm::Handle<edm::View<pat::Muon>> muons_;
    event.getByLabel("slimmedMuons",muons_);
    
    
    
    bool skipEM_                    = true;
    double skipEMfractionThreshold_ = 0.9;
    bool skipMuons_                 = true;
    
    double jetCorrEtaMax_           = 9.9;
    double type1JetPtThreshold_     = 10.0;
    
    double corrEx    = 0;
    double corrEy    = 0;
    double corrSumEt = 0;
    
    for (const pat::Jet &jet : *jets_) {
        
        double emEnergyFraction = jet.chargedEmEnergyFraction() + jet.neutralEmEnergyFraction();
        if ( skipEM_ && emEnergyFraction > skipEMfractionThreshold_ ) continue;
        
        reco::Candidate::LorentzVector rawJetP4 = jet.correctedP4(0);
        double corr = getJEC(rawJetP4, jet, jetCorrEtaMax_, jetCorrLabel_);
        
        if ( skipMuons_ && jet.muonMultiplicity() != 0 ) {
            
            for (const pat::Muon &muon : *muons_) {
                if( !muon.isGlobalMuon() && !muon.isStandAloneMuon() ) continue;
                TLorentzVector muonV; muonV.SetPtEtaPhiE(muon.p4().pt(),muon.p4().eta(),muon.p4().phi(),muon.p4().e());
                TLorentzVector jetV; jetV.SetPtEtaPhiE(jet.p4().pt(),jet.p4().eta(),jet.p4().phi(),jet.p4().e());
                if( muonV.DeltaR(jetV) < 0.5 ){
                    reco::Candidate::LorentzVector muonP4 = muon.p4();
                    rawJetP4 -= muonP4;
                }
            }
        }
        
        reco::Candidate::LorentzVector corrJetP4 = corr*rawJetP4;
        
        if ( corrJetP4.pt() > type1JetPtThreshold_ ) {
            reco::Candidate::LorentzVector tmpP4 = jet.correctedP4(0);
            corr = getJECOffset(tmpP4, jet, jetCorrEtaMax_, offsetCorrLabel_);
            reco::Candidate::LorentzVector rawJetP4offsetCorr = corr*rawJetP4;
            
            corrEx    -= (corrJetP4.px() - rawJetP4offsetCorr.px());
            corrEy    -= (corrJetP4.py() - rawJetP4offsetCorr.py());
            corrSumEt += (corrJetP4.Et() - rawJetP4offsetCorr.Et());
        }
    }
    TypeICorrMap_["corrEx"]    = corrEx;
    TypeICorrMap_["corrEy"]    = corrEy;
    TypeICorrMap_["corrSumEt"] = corrSumEt;
}

*/
bool PKUTreeMaker::hasMatchedPromptElectron(const reco::SuperClusterRef &sc, const edm::Handle<edm::View<pat::Electron> > &eleCol, const edm::Handle<reco::ConversionCollection> &convCol, const math::XYZPoint &beamspot,  float lxyMin, float probMin, unsigned int nHitsBeforeVtxMax) {
    //check if a given SuperCluster matches to at least one GsfElectron having zero expected inner hits
    //and not matching any conversion in the collection passing the quality cuts
    if (sc.isNull()) return false;
    for (edm::View<pat::Electron>::const_iterator it = eleCol->begin(); it!=eleCol->end(); ++it) {
        //match electron to supercluster
        if (it->superCluster()!=sc) continue;
        //check expected inner hits
        if (it->gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) > 0) continue;
        //check if electron is matching to a conversion
//        if (ConversionTools::hasMatchedConversion(*it,convCol,beamspot)) continue;
        return true;
    }
    return false;
}


PKUTreeMaker::~PKUTreeMaker()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}

//-------------------------------------------------------------------------------------------------------------------------------------//
//
// member functions
//

// ------------ method called for each event  ------------
void
PKUTreeMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   nevent = iEvent.eventAuxiliary().event();
   run    = iEvent.eventAuxiliary().run();
   ls     = iEvent.eventAuxiliary().luminosityBlock();

   edm::Handle<edm::View<reco::Candidate> > leptonicVs;
   iEvent.getByLabel(leptonicVSrc_.c_str(), leptonicVs);
    

    iEvent.getByToken(rhoToken_      , rho_     );
    double fastJetRho = *(rho_.product());
    useless = fastJetRho;
  
   edm::Handle<edm::View<pat::Jet> > ak4jets;
   iEvent.getByLabel(ak4jetsSrc_.c_str(), ak4jets);

   edm::Handle<edm::View<pat::Photon> > photons;
   iEvent.getByLabel(photonSrc_.c_str(), photons);

   edm::Handle<edm::View<reco::Candidate> > metHandle;
   iEvent.getByLabel(metSrc_.c_str(), metHandle);
  
   edm::Handle<edm::View<reco::GenParticle> > genParticles;//define genParticle
   iEvent.getByLabel(InputTag("prunedGenParticles"), genParticles);

   edm::Handle<edm::View<pat::Muon>> mus;
   iEvent.getByLabel("slimmedMuons",mus);

   edm::Handle<edm::View<pat::Electron>> eles;
   iEvent.getByLabel("slimmedElectrons",eles);


 
   //StringCutObjectSelector<T>
   
   setDummyValues(); //Initalize variables with dummy values
   

//-------------------------------------------------------------------------------------------------------------------------------------//
            //*****************************************************************//
            //************************* Gen Level Information******************//
/*
	  for(size_t ik=0; ik<genParticles->size();ik++)
	{
		if( (*genParticles)[ik].pdgId()==5100039 ) // && (*genParticles)[ik].status()==3)//graviton
			}//end of graviton daughter loop
		}//end of graviton
         }
           
*/
//-------------------------------------------------------------------------------------------------------------------------------------//

           

       const reco::Candidate& leptonicV = leptonicVs->at(0);
       const reco::Candidate& metCand = metHandle->at(0);
       
           edm::Handle<reco::VertexCollection> vertices;
           iEvent.getByLabel("offlineSlimmedPrimaryVertices", vertices);
           if (vertices->empty()) return; // skip the event if no PV found
	   nVtx = vertices->size();
           reco::VertexCollection::const_iterator firstGoodVertex = vertices->end();
           for (reco::VertexCollection::const_iterator vtx = vertices->begin(); vtx != vertices->end(); ++vtx) {
               // Replace isFake() for miniAOD because it requires tracks and miniAOD vertices don't have tracks:
               // Vertex.h: bool isFake() const {return (chi2_==0 && ndof_==0 && tracks_.empty());}
               if (  /*!vtx->isFake() &&*/ 
                     !(vtx->chi2()==0 && vtx->ndof()==0) 
	             &&  vtx->ndof()>=4. && vtx->position().Rho()<=2.0
	             && fabs(vtx->position().Z())<=24.0) {
                  firstGoodVertex = vtx;
                  break;
               }           
           }
           if ( firstGoodVertex==vertices->end() ) return; // skip event if there are no good PVs

    //*****************************************************************//
    //************************* MET **********************//
    //*****************************************************************//
 /*   edm::Handle<pat::METCollection>  METs_;
    bool defaultMET = iEvent.getByToken(metInputToken_ , METs_ );
    edm::Handle<pat::METCollection>  reclusteredMETs_;
    bool reclusteredMET = iEvent.getByToken(reclusteredmetInputToken_ , reclusteredMETs_ );
    edm::Handle<edm::View<reco::PFMET> >     pfMET_ ;
    if ( METsRawLabel_.label() != "" ) iEvent.getByLabel(METsRawLabel_ , pfMET_ );
    if(defaultMET){
        
        addTypeICorr(iEvent);
        for (const pat::MET &met : *METs_) {

            const float rawPt  = met.shiftedPt(pat::MET::NoShift, pat::MET::Raw);
            const float rawPhi = met.shiftedPhi(pat::MET::NoShift, pat::MET::Raw);
            const float rawSumEt = met.shiftedSumEt(pat::MET::NoShift, pat::MET::Raw);

            TVector2 rawMET_;
            rawMET_.SetMagPhi (rawPt, rawPhi );
            
            
            Double_t rawPx = rawMET_.Px();
            Double_t rawPy = rawMET_.Py();
            Double_t rawEt = std::hypot(rawPx,rawPy);
            
            METraw_et = rawEt;
            METraw_phi = rawPhi;
            METraw_sumEt = rawSumEt;
            
            double pxcorr = rawPx+TypeICorrMap_["corrEx"];
            double pycorr = rawPy+TypeICorrMap_["corrEy"];
            double et     = std::hypot(pxcorr,pycorr);
            
            double sumEtcorr = rawSumEt+TypeICorrMap_["corrSumEt"];
            TLorentzVector corrmet; corrmet.SetPxPyPzE(pxcorr,pycorr,0.,et);
            useless = sumEtcorr;
            useless = rawEt;
            MET_et = et;
            MET_phi = corrmet.Phi();
            MET_sumEt = sumEtcorr;
            MET_corrPx = TypeICorrMap_["corrEx"];
            MET_corrPy = TypeICorrMap_["corrEy"];
        }
    }
    
    else if (reclusteredMET){
        METraw_et = (pfMET_->front() ).et();
        METraw_phi = (pfMET_->front() ).phi();
        METraw_sumEt = (pfMET_->front() ).sumEt();
        
        for (const pat::MET &met : *reclusteredMETs_) {
            double MET_et = met.et();
            MET_et = met.et();
            MET_phi = met.phi();
            MET_sumEt = met.sumEt();
            useless = MET_et;
        }
    }
*/
  

          //*****************************************************************//
           //************************* ID for electrons **********************//
           //*****************************************************************//
                 if( leptonicV.daughter(0)->isElectron()||leptonicV.daughter(1)->isElectron() ) {
                       const pat::Electron *el1 = leptonicV.daughter(0)->isElectron() ? 
                                                  (pat::Electron*)leptonicV.daughter(0):
                                                  (pat::Electron*)leptonicV.daughter(1);
                    if (el1->gsfTrack().isNonnull()){
                        reco::GsfElectron::PflowIsolationVariables pfIso1 = el1->pfIsolationVariables();
                        ptel1          = el1->pt();
                        etaSC1         = el1->superCluster()->eta();
                        dEtaIn1        = el1->deltaEtaSuperClusterTrackAtVtx();
                        dPhiIn1        = el1->deltaPhiSuperClusterTrackAtVtx();
                        hOverE1        = el1->hcalOverEcal();
                        full5x5_sigma1 = el1->full5x5_sigmaIetaIeta();
                        ooEmooP1       = el1->ecalEnergy() && std::isfinite(el1->ecalEnergy()) ? 
                                         fabs(1.0/el1->ecalEnergy() - el1->eSuperClusterOverP()/el1->ecalEnergy() ) : 1e9;
                        double absiso1 = pfIso1.sumChargedHadronPt + std::max(0.0, pfIso1.sumNeutralHadronEt + pfIso1.sumPhotonEt - 0.5*pfIso1.sumPUPt );
                        relIso1        = absiso1/el1->pt();
                        d01            = (-1)*el1->gsfTrack()->dxy(firstGoodVertex->position());   
                        dz1            = el1->gsfTrack()->dz(firstGoodVertex->position());
                        missingHits1   = el1->gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS);
            	        passConVeto1   = el1->passConversionVeto();
                        el1passID      = ( el1->electronID(electronIdTag_.encode())>0.5 );

                        if( fabs(ooEmooP1)<0.05 && fabs(d01)<0.02 && fabs(dz1)<0.1 && missingHits1==0 && passConVeto1==1 && relIso1<0.10 && ( ( fabs(etaSC1)<=1.479 && full5x5_sigma1<0.01 && fabs(dPhiIn1)<0.03 && fabs(dEtaIn1)<0.004 && hOverE1<0.12 ) ||   (fabs(etaSC1)>1.479 && fabs(etaSC1)<2.5 && full5x5_sigma1<0.03 && fabs(dPhiIn1)<0.02 && fabs(dEtaIn1)<0.005 && ooP1)<0.05 && fabs(d01)<0.02 && fabs(dz1)<0.1 && missingHits1<2 && passConVeto1==1 && relIso1<0.15 && ( ( fabs(etaSC1)<=1.479 && full5x5_sigma1<0.01 && fabs(dPhiIn1)<0.06 && fabs(dEtaIn1)<0.004 && hOverE1<0.12 ) ||   (fabs(etaSC1)>1.479 && fabs(etaSC1)<2.5 && full5x5_sigma1<0.03 && fabs(dPhiIn1)<0.03 && fabs(dEtaIn1)<0.007 && hOverE1<0.10 )) ) {
                        el1mediumID =1 ; 
                        }

                        if( fabs(ooEmooP1)<0.05 && fabs(d01)<0.02 && fabs(dz1)<0.1 && missingHits1<2 && passConVeto1==1 && relIso1<0.15 && ( ( fabs(etaSC1)<=1.479 && full5x5_sigma1<0.01 && fabs(dPhiIn1)<0.15 && fabs(dEtaIn1)<0.007 && hOverE1<0.12 ) ||   (fabs(etaSC1)>1.479 && fabs(etaSC1)<2.5 && full5x5_sigma1<0.03 && fabs(dPhiIn1)<0.10 && fabs(dEtaIn1)<0.009 && hOverE1<0.10 )) ) {
                        el1looseID =1 ; 
                        }

}}
//---------------------------



//                 if( abs(leptonicV.daughter(0)->pdgId())==13 || abs(leptonicV.daughter(1)->pdgId())==13 ) {
                   if( leptonicV.daughter(0)->isMuon()||leptonicV.daughter(1)->isMuon()){

                       const pat::Muon *mu1 = abs(leptonicV.daughter(0)->pdgId())==13 ? 
                                                  (pat::Muon*)leptonicV.daughter(0):
                                                  (pat::Muon*)leptonicV.daughter(1);

                       double absiso1_mu = (mu1->pfIsolationR04().sumChargedHadronPt+ std::max(0.0,mu1->pfIsolationR04().sumNeutralHadronEt+mu1->pfIsolationR04().sumPhotonEt-0.5*mu1->pfIsolationR04().sumPUPt))/mu1->pt();                        


                       if(mu1-> isGlobalMuon() && mu1->isPFMuon() && (mu1->globalTrack()->normalizedChi2())<10 && (mu1->globalTrack()->hitPattern().numberOfValidMuonHits())>0 && (mu1->numberOfMatchedStations())>1 && fabs(mu1->muonBestTrack()->dxy(firstGoodVertex->position()))<0.2 && fabs(mu1->muonBestTrack()->dz(firstGoodVertex->position()))<0.5 && (mu1->numberOfMatchedStations())>1 && (mu1->globalTrack()->hitPattern().trackerLayersWithMeasurement())>5 && fabs(absiso1_mu)<0.12) {
                        mu1tightID =1 ;
                        } 


                  }            

//---------------------------
                


       /// For the time being, set these to 1
       triggerWeight=1.0;
       pileupWeight=1.0;

       double targetEvents = targetLumiInvPb_*crossSectionPb_;
       lumiWeight = targetEvents/originalNEvents_;

       ptVlep       = leptonicV.pt();
       yVlep        = leptonicV.eta();
       phiVlep      = leptonicV.phi();
       massVlep     = leptonicV.mass();
       mtVlep       = leptonicV.mt();
   

       ptlep1       = leptonicV.daughter(1)->pt();
       etalep1      = leptonicV.daughter(1)->eta();
       philep1      = leptonicV.daughter(1)->phi();
       if(leptonicV.daughter(0)->isElectron()||leptonicV.daughter(0)->isMuon() ) {
       ptlep1       = leptonicV.daughter(0)->pt();
       etalep1      = leptonicV.daughter(0)->eta();
       philep1      = leptonicV.daughter(0)->phi(); }


       met          = metCand.pt();
       metPhi       = metCand.phi();

//       deltaRleplep = deltaR(etalep1,philep1,etalep2,philep2);
//       deltaRlepjet = std::min(drl1j,drl2j);
//       delPhilepmet = deltaPhi(philep1, metPhi);

        lep          = std::max(abs(leptonicV.daughter(0)->pdgId()), abs(leptonicV.daughter(1)->pdgId()));


            //************************* Photon Jets Information******************//

         double rhoVal_;
         rhoVal_=-99.;
         //edm::Handle<double> rho;
        // iEvent.getByLabel("fixedGridRhoFastjetAll",rho);
         rhoVal_ = *rho_;


         photonet=-100.;  iphoton=-1;

            for (size_t ip=0; ip<photons->size();ip++)
         {
            if(ip<6)  {
                photon_pt[ip] = (*photons)[ip].pt();
                photon_eta[ip] = (*photons)[ip].eta();
                photon_phi[ip] = (*photons)[ip].phi();
                photon_e[ip] = (*photons)[ip].energy();
            }

            int istightphoton=0;

            // photon id
            // ??ele convert veto ?? (*photons)[ip].superCluster().eta()
             edm::Handle<edm::View<pat::Electron> > electrons;
             iEvent.getByToken(electronToken_, electrons);
             
             edm::Handle<reco::BeamSpot> beamSpot;
             iEvent.getByToken(beamSpotToken_,beamSpot);
             
             edm::Handle<std::vector<reco::Conversion> > conversions;
             iEvent.getByToken(conversionsToken_,conversions);
             
             edm::Handle<edm::View<pat::Photon> > photons;
             iEvent.getByToken(photonToken_, photons);
             
             bool passEleVeto = (!hasMatchedPromptElectron((*photons)[ip].superCluster(),electrons, conversions, beamSpot->position() ) );
             std::cout << "PhoPassEleVeto " << passEleVeto << std::endl;

            double phoiso=std::max((*photons)[ip].photonIso()-rhoVal_*EApho(fabs((*photons)[ip].eta())),0.0);
            double chiso=std::max((*photons)[ip].chargedHadronIso()-rhoVal_*EAch(fabs((*photons)[ip].eta())),0.0);
            double nhiso=std::max((*photons)[ip].neutralHadronIso()-rhoVal_*EAnh(fabs((*photons)[ip].eta())),0.0);

            if((*photons)[ip].isEB() && (*photons)[ip].hadTowOverEm()<0.011 && (*photons)[ip].sigmaIetaIeta()<0.0099 && chiso<1.86 && nhiso<(2.64 + 0.0025*(*photons)[ip].pt()) && phoiso<(1.20+0.001*(*photons)[ip].pt())) {istightphoton=1;}
            if((*photons)[ip].isEE() && (*photons)[ip].hadTowOverEm()<0.015 && (*photons)[ip].sigmaIetaIeta()<0.0263 && chiso<1.68 && nhiso<(4.42 + 0.0118*(*photons)[ip].pt()) && phoiso<(1.03+0.0059*(*photons)[ip].pt())) {istightphoton=1;}
            //std::cout<<"tight photon "<<istightphoton<<std::endl;

            //
             if(istightphoton==1 && deltaR(photon_eta[ip],photon_phi[ip],etalep1,philep1) > 0.5) {
                 if(photon_pt[ip]>photonet) {iphoton=ip;}
             }     
         }

         if(iphoton>-1) {
               photonet=(*photons)[iphoton].pt();
               photoneta=(*photons)[iphoton].eta();
               photonphi=(*photons)[iphoton].phi();
               photone=(*photons)[iphoton].energy();
               drla=deltaR(photon_eta[iphoton],photon_phi[iphoton],etalep1,philep1);
         }  


//std::cout<<iphoton<<" "<<photonet<<std::endl;

            //************************* AK4 Jets Information******************//
           
//         for ( auto jet =ak4jets->begin(); jet != ak4jets->end(); ++jet) 

    Int_t jetindexphoton12[2] = {-1,-1};
    Int_t jetindexphoton12_jesup[2] = {-1,-1};
    Int_t jetindexphoton12_jesdown[2] = {-1,-1};
//    Int_t jetindexphoton12_jerup[2] = {0,1};
  //  Int_t jetindexphoton12_jerdown[2] = {-1,-1};


         std::vector<JetCorrectorParameters> vPar;
//         vPar.clear();
         for ( std::vector<std::string>::const_iterator payloadBegin = jecAK4Labels_.begin(), payloadEnd = jecAK4Labels_.end(), ipayload = payloadBegin; ipayload != payloadEnd; ++ipayload ) {
         JetCorrectorParameters pars(*ipayload);
         vPar.push_back(pars);
        }

        jecAK4_ = new FactorizedJetCorrector(vPar);


        int nujets=0 ;
        double tmpjetptcut=20.0;
        std::vector<TLorentzVector*> jets;
        std::vector<TLorentzVector*> ak4jet_jesup_p4;
    std::vector<TLorentzVector*> ak4jet_jesdown_p4;
    std::vector<TLorentzVector*> ak4jet_jerup_p4;
    std::vector<TLorentzVector*> ak4jet_jerdown_p4;
    
//################Jet Correction##########################

            for (size_t ik=0; ik<ak4jets->size();ik++)
         {

            reco::Candidate::LorentzVector uncorrJet = (*ak4jets)[ik].correctedP4(0);
            jecAK4_->setJetEta( uncorrJet.eta() );
            jecAK4_->setJetPt ( uncorrJet.pt() );
            jecAK4_->setJetE ( uncorrJet.energy() );
            jecAK4_->setRho ( rhoVal_ );
            jecAK4_->setNPV ( vertices->size() );
            jecAK4_->setJetA ( (*ak4jets)[ik].jetArea() );
            double corr = jecAK4_->getCorrection();

            if(corr*uncorrJet.pt()>tmpjetptcut) {
            TLorentzVector *dummy = new TLorentzVector(0,0,0,0);    
            dummy->SetPtEtaPhiE(corr*uncorrJet.pt(), uncorrJet.eta(), uncorrJet.phi(), corr*uncorrJet.energy());
            jets.push_back(dummy);
            ++nujets;
            }   

            if(ik<6)  {   
                ak4jet_pt[ik] =  corr*uncorrJet.pt();
                ak4jet_eta[ik] = (*ak4jets)[ik].eta();
                ak4jet_phi[ik] = (*ak4jets)[ik].phi();
                ak4jet_e[ik] =   corr*uncorrJet.energy();
                ak4jet_csv[ik] = (*ak4jets)[ik].bDiscriminator("pfCombinedSecondaryVertexBJetTags");
                ak4jet_icsv[ik] = (*ak4jets)[ik].bDiscriminator("combinedInclusiveSecondaryVertexV2BJetTags");
            }
          }

//################JES Down and Up##########################

    for (size_t ik=0; ik<ak4jets->size();ik++)
    {
    jesAK4_ = new getJesUnc();
jesAK4_->setJetPtEta(ak4jet_pt[ik], ak4jet_eta[ik]);
        double jesuncert = jesAK4_->getjesunc();
        
        double SF_ak4jet_jesup = 1.0 + jesuncert;
        double SF_ak4jet_jesdown = 1.0 - jesuncert;
        
        
        ak4jet_jesup_pt_[ik] = SF_ak4jet_jesup * ak4jet_pt[ik];
        ak4jet_jesdown_pt_[ik] = SF_ak4jet_jesdown * ak4jet_pt[ik];
        ak4jet_jesup_e_[ik] = SF_ak4jet_jesup * ak4jet_e[ik];
        ak4jet_jesdown_e_[ik] = SF_ak4jet_jesdown * ak4jet_e[ik];
        
       // if(SF_ak4jet_jesup * ak4jet_pt[ik]>tmpjetptcut){
        TLorentzVector *dummy = new TLorentzVector(0,0,0,0);
        dummy->SetPtEtaPhiE(ak4jet_jesup_pt_[ik], ak4jet_eta[ik], ak4jet_phi[ik], ak4jet_jesup_e_[ik]);
        ak4jet_jesup_p4.push_back(dummy);

      //  }
      //  if(SF_ak4jet_jesdown * ak4jet_pt[ik]>tmpjetptcut){
        dummy->SetPtEtaPhiE(ak4jet_jesdown_pt_[ik], ak4jet_eta[ik], ak4jet_phi[ik], ak4jet_jesdown_e_[ik]);
        ak4jet_jesdown_p4.push_back(dummy);
     //   }
    }
//################JER Down and Up##########################
    
    for (size_t ik=0; ik<ak4jets->size();ik++)
    {
        jetresolutiontable_ = new LOTable();
        jetresolutiontable_->LoadTable("TagJetResolution.txt");
        double jeruncert = jetresolutiontable_->GetValue(ak4jet_pt[ik],fabs(ak4jet_eta[ik]));
        //std::cout<<jeruncert<<std::endl;
        jerAK4_ = new jetResolution(jeruncert);
        jerAK4_->SetSeed(ak4jet_phi[ik]);

        jerAK4_->SetJetPtEta(ak4jet_pt[ik],fabs(ak4jet_eta[ik]));
        ak4jet_jerup_pt_[ik] = jerAK4_->GetModifiedJetPt(0);
        ak4jet_jerdown_pt_[ik] = jerAK4_->GetModifiedJetPt(1);
        ak4jet_pt[ik] = jerAK4_->GetModifiedJetPt(2);
        std::cout<<ak4jet_jerup_pt_[ik]<<" "<<ak4jet_jerdown_pt_[ik]<<" "<<ak4jet_pt[ik]<<std::endl;
        
        double SF_ak4jet_jerup = ak4jet_jerup_pt_[ik] / ak4jet_pt[ik];
        double SF_ak4jet_jerdown = ak4jet_jerdown_pt_[ik] / ak4jet_pt[ik];
        
        
        ak4jet_jerup_e_[ik] = SF_ak4jet_jerup * ak4jet_e[ik];
        ak4jet_jerdown_e_[ik] = SF_ak4jet_jerdown * ak4jet_e[ik];
        
        //  if(ak4jet_jerup_pt_[ik]>tmpjetptcut){
        TLorentzVector *dummy = new TLorentzVector(0,0,0,0);
        dummy->SetPtEtaPhiE(ak4jet_jerup_pt_[ik], ak4jet_eta[ik], ak4jet_phi[ik], ak4jet_jerup_e_[ik]);
        ak4jet_jerup_p4.push_back(dummy);
   // }
       // if(ak4jet_jerdown_pt_[ik]>tmpjetptcut){
        dummy->SetPtEtaPhiE(ak4jet_jerdown_pt_[ik], ak4jet_eta[ik], ak4jet_phi[ik], ak4jet_jerdown_e_[ik]);
        ak4jet_jerdown_p4.push_back(dummy);
       // }
            }


  
    
    //***********************normal jet****************************//

    sort (jets.begin (), jets.end (), mysortPt);

           for (size_t i=0;i<jets.size();i++) {
             if(iphoton>0) {
               double drtmp1=deltaR(jets.at(i)->Eta(), jets.at(i)->Phi(), photoneta,photonphi);

              if(drtmp1>0.5 && jetindexphoton12[0]==-1&&jetindexphoton12[1]==-1) {
                     jetindexphoton12[0] = i;
                     continue;
              }
              if(drtmp1>0.5 && jetindexphoton12[0]!=-1&&jetindexphoton12[1]==-1) {
                     jetindexphoton12[1] = i;
                     continue;
              }
            }
         }


         if(jetindexphoton12[0]>-1 && jetindexphoton12[1]>-1) {
            jet1pt=jets[jetindexphoton12[0]]->Pt();
            jet1eta=jets[jetindexphoton12[0]]->Eta();
            jet1phi=jets[jetindexphoton12[0]]->Phi();
            jet1e=jets[jetindexphoton12[0]]->E();
//            jet1csv=(*ak4jets)[jetindexphoton12[0]].bDiscriminator("pfCombinedSecondaryVertexBJetTags");
//            jet1icsv=(*ak4jets)[jetindexphoton12[0]].bDiscriminator("combinedInclusiveSecondaryVertexV2BJetTags");
            jet2pt=jets[jetindexphoton12[1]]->Pt();
            jet2eta=jets[jetindexphoton12[1]]->Eta();
            jet2phi=jets[jetindexphoton12[1]]->Phi();
            jet2e=jets[jetindexphoton12[1]]->E();
//            jet2csv=(*ak4jets)[jetindexphoton12[1]].bDiscriminator("pfCombinedSecondaryVertexBJetTags");
//            jet2icsv=(*ak4jets)[jetindexphoton12[1]].bDiscriminator("combinedInclusiveSecondaryVertexV2BJetTags");

            drj1a=deltaR(jet1eta,jet1phi,photoneta,photonphi);
            drj2a=deltaR(jet2eta,jet2phi,photoneta,photonphi);
            drj1l=deltaR(jet1eta,jet1phi,etalep1,philep1);
            drj2l=deltaR(jet2eta,jet2phi,etalep1,philep1);
            
            TLorentzVector j1p4;
            j1p4.SetPtEtaPhiE(jet1pt, jet1eta, jet1phi, jet1e);
            TLorentzVector j2p4;
            j2p4.SetPtEtaPhiE(jet2pt, jet2eta, jet2phi, jet2e);
            TLorentzVector photonp4;
            photonp4.SetPtEtaPhiE(photonet, photoneta, photonphi, photone);
            TLorentzVector vp4;
            vp4.SetPtEtaPhiE(leptonicV.pt(), leptonicV.eta(), leptonicV.phi(), leptonicV.energy());

            j1metPhi=fabs(jet1phi-metPhi);
            if(j1metPhi>Pi) {j1metPhi=2.0*Pi-j1metPhi;}
            

            j2metPhi=fabs(jet2phi-metPhi);
            if(j2metPhi>Pi) {j2metPhi=2.0*Pi-j2metPhi;}

 
            Mjj=(j1p4 + j2p4).M();
            deltaeta = fabs(jet1eta - jet2eta);
            zepp = fabs((vp4+photonp4).Rapidity() - (j1p4.Rapidity() + j2p4.Rapidity())/ 2.0); 
         }
    //***********************jes up jet****************************//
    sort (ak4jet_jesup_p4.begin (), ak4jet_jesup_p4.end (), mysortPt);
    
    for (size_t i=0;i<ak4jet_jesup_p4.size();i++) {
        if(iphoton>0) {
            double drtmp1=deltaR(ak4jet_jesup_p4.at(i)->Eta(), ak4jet_jesup_p4.at(i)->Phi(), photoneta,photonphi);
            
            if(drtmp1>0.5 && jetindexphoton12_jesup[0]==-1&&jetindexphoton12_jesup[1]==-1) {
                jetindexphoton12_jesup[0] = i;
                continue;
            }
            if(drtmp1>0.5 && jetindexphoton12_jesup[0]!=-1&&jetindexphoton12_jesup[1]==-1) {
                jetindexphoton12_jesup[1] = i;
                continue;
            }
        }
    }
    
    
    if(jetindexphoton12_jesup[0]>-1 && jetindexphoton12_jesup[1]>-1) {
        jet1_jesup_pt=ak4jet_jesup_p4[jetindexphoton12_jesup[0]]->Pt();
        jet1_jesup_eta=ak4jet_jesup_p4[jetindexphoton12_jesup[0]]->Eta();
        jet1_jesup_phi=ak4jet_jesup_p4[jetindexphoton12_jesup[0]]->Phi();
        jet1_jesup_e=ak4jet_jesup_p4[jetindexphoton12_jesup[0]]->E();
        
        jet2_jesup_pt=ak4jet_jesup_p4[jetindexphoton12_jesup[1]]->Pt();
        jet2_jesup_eta=ak4jet_jesup_p4[jetindexphoton12_jesup[1]]->Eta();
        jet2_jesup_phi=ak4jet_jesup_p4[jetindexphoton12_jesup[1]]->Phi();
        jet2_jesup_e=ak4jet_jesup_p4[jetindexphoton12_jesup[1]]->E();
        
        TLorentzVector j1_jesup_p4;
        j1_jesup_p4.SetPtEtaPhiE(jet1_jesup_pt, jet1_jesup_eta, jet1_jesup_phi, jet1_jesup_e);
        TLorentzVector j2_jesup_p4;
        j2_jesup_p4.SetPtEtaPhiE(jet2_jesup_pt, jet2_jesup_eta, jet2_jesup_phi, jet2_jesup_e);
        
        jet_jesup_mass = (j1_jesup_p4 + j2_jesup_p4).M();
    }

    //***********************jes down jet****************************//
    
    
    sort (ak4jet_jesdown_p4.begin (), ak4jet_jesdown_p4.end (), mysortPt);
   size_t num=ak4jet_jesdown_p4.size();
    std::cout<<" num "<<num<<std::endl; 
    for (size_t i=0;i<ak4jet_jesdown_p4.size();i++) {
    std::cout<<" iphoton "<<iphoton<<std::endl;       
 if(iphoton>0) {
            double drtmp1=deltaR(ak4jet_jesdown_p4.at(i)->Eta(), ak4jet_jesdown_p4.at(i)->Phi(), photoneta,photonphi);
                std::cout<<" drtmp1 "<<drtmp1<<std::endl;
            if(drtmp1>0.5 && jetindexphoton12_jesdown[0]==-1&&jetindexphoton12_jesdown[1]==-1) {
                jetindexphoton12_jesdown[0] = i;
                continue;
            }
            if(drtmp1>0.5 && jetindexphoton12_jesdown[0]!=-1&&jetindexphoton12_jesdown[1]==-1) {
                jetindexphoton12_jesdown[1] = i;
                continue;
            }
        }
    }
    
    
    if(jetindexphoton12_jesdown[0]>-1 && jetindexphoton12_jesdown[1]>-1) {
        jet1_jesdown_pt=ak4jet_jesdown_p4[jetindexphoton12_jesdown[0]]->Pt();
        jet1_jesdown_eta=ak4jet_jesdown_p4[jetindexphoton12_jesdown[0]]->Eta();
        jet1_jesdown_phi=ak4jet_jesdown_p4[jetindexphoton12_jesdown[0]]->Phi();
        jet1_jesdown_e=ak4jet_jesdown_p4[jetindexphoton12_jesdown[0]]->E();
std::cout<<jet1_jesdown_pt<<" "<<jet1_jesdown_eta<<" "<<jet1_jesdown_phi<<" "<<jet1_jesdown_e<<std::endl;        
        jet2_jesdown_pt=ak4jet_jesdown_p4[jetindexphoton12_jesdown[1]]->Pt();
        jet2_jesdown_eta=ak4jet_jesdown_p4[jetindexphoton12_jesdown[1]]->Eta();
        jet2_jesdown_phi=ak4jet_jesdown_p4[jetindexphoton12_jesdown[1]]->Phi();
        jet2_jesdown_e=ak4jet_jesdown_p4[jetindexphoton12_jesdown[1]]->E();
        std::cout<<jet2_jesdown_pt<<" "<<jet2_jesdown_eta<<" "<<jet2_jesdown_phi<<" "<<jet2_jesdown_e<<std::endl;
        TLorentzVector j1_jesdown_p4;
        j1_jesdown_p4.SetPtEtaPhiE(jet1_jesdown_pt, jet1_jesdown_eta, jet1_jesdown_phi, jet1_jesdown_e);
        TLorentzVector j2_jesdown_p4;
        j2_jesdown_p4.SetPtEtaPhiE(jet2_jesdown_pt, jet2_jesdown_eta, jet2_jesdown_phi, jet2_jesdown_e);
        
        jet_jesdown_mass = (j1_jesdown_p4 + j2_jesdown_p4).M();
    }

    //***********************jer up jet****************************//
    sort (ak4jet_jerup_p4.begin (), ak4jet_jerup_p4.end (), mysortPt);
    
  /*  size_t num=ak4jet_jerup_p4.size();
    std::cout<<" num "<<num<<std::endl;
    
    for (size_t i=0;i<ak4jet_jerup_p4.size();i++) {
    std::cout<<" iphoton "<<iphoton<<std::endl;        
    if(iphoton>0) {
            double drtmp1=deltaR(ak4jet_jerup_p4.at(i)->Eta(), ak4jet_jerup_p4.at(i)->Phi(), photoneta,photonphi);
       std::cout<<" drtmp1 "<<drtmp1<<std::endl;            
            if(drtmp1>0.5 && jetindexphoton12_jerup[0]==-1&&jetindexphoton12_jerup[1]==-1) {
                jetindexphoton12_jerup[0] = i;
                continue;
            }
            if(drtmp1>0.5 && jetindexphoton12_jerup[0]!=-1&&jetindexphoton12_jerup[1]==-1) {
                jetindexphoton12_jerup[1] = i;
                continue;
           }
}
 std::cout<<"333"<<std::endl;    
}
    std::cout<<"444"<<std::endl;
  */
        if(num1>1){
        if(jetindexphoton12[0]>-1 && jetindexphoton12[1]>-1) {
        jet1_jerup_pt=ak4jet_jerup_p4[jetindexphoton12[0]]->Pt();
        jet1_jerup_eta=ak4jet_jerup_p4[jetindexphoton12[0]]->Eta();
        jet1_jerup_phi=ak4jet_jerup_p4[jetindexphoton12[0]]->Phi();
        jet1_jerup_e=ak4jet_jerup_p4[jetindexphoton12[0]]->E();
        
        jet2_jerup_pt=ak4jet_jerup_p4[jetindexphoton12[1]]->Pt();
        jet2_jerup_eta=ak4jet_jerup_p4[jetindexphoton12[1]]->Eta();
        jet2_jerup_phi=ak4jet_jerup_p4[jetindexphoton12[1]]->Phi();
        jet2_jerup_e=ak4jet_jerup_p4[jetindexphoton12[1]]->E();
        
        TLorentzVector j1_jerup_p4;
        j1_jerup_p4.SetPtEtaPhiE(jet1_jerup_pt, jet1_jerup_eta, jet1_jerup_phi, jet1_jerup_e);
        TLorentzVector j2_jerup_p4;
        j2_jerup_p4.SetPtEtaPhiE(jet2_jerup_pt, jet2_jerup_eta, jet2_jerup_phi, jet2_jerup_e);
        
        jet_jerup_mass = (j1_jerup_p4 + j2_jerup_p4).M();
    }
        }
    //***********************jer down jet****************************
    
    
    sort (ak4jet_jerdown_p4.begin (), ak4jet_jerdown_p4.end (), mysortPt);
   size_t num1=ak4jet_jerdown_p4.size();
    std::cout<<" num1 "<<num1<<std::endl; 
/*    for (size_t i=0;i<ak4jet_jerdown_p4.size();i++) {
std::cout<<" iphoton "<<iphoton<<std::endl;        
if(iphoton>0) {
            double drtmp1=deltaR(ak4jet_jerdown_p4.at(i)->Eta(), ak4jet_jerdown_p4.at(i)->Phi(), photoneta,photonphi);
           std::cout<<" drtmp1 "<<drtmp1<<std::endl; 
            if(drtmp1>0.5 && jetindexphoton12_jerdown[0]==-1&&jetindexphoton12_jerdown[1]==-1) {
                jetindexphoton12_jerdown[0] = i;
                continue;
            }
            if(drtmp1>0.5 && jetindexphoton12_jerdown[0]!=-1&&jetindexphoton12_jerdown[1]==-1) {
                jetindexphoton12_jerdown[1] = i;
                continue;
            }
        }
    }
    
  */
    if(num1>1){
        std::cout<<jetindexphoton12[0]<<" "<<jetindexphoton12[1]<<std::endl;
    if(jetindexphoton12[0]>-1 && jetindexphoton12[1]>-1) {
        jet1_jerdown_pt=ak4jet_jerdown_p4[jetindexphoton12[0]]->Pt();
        jet1_jerdown_eta=ak4jet_jerdown_p4[jetindexphoton12[0]]->Eta();
        jet1_jerdown_phi=ak4jet_jerdown_p4[jetindexphoton12[0]]->Phi();
        jet1_jerdown_e=ak4jet_jerdown_p4[jetindexphoton12[0]]->E();
       std::cout<<jet1_jerdown_pt<<" "<<jet1_jerdown_eta<<" "<<jet1_jerdown_phi<<" "<<jet1_jerdown_e<<std::endl;
        jet2_jerdown_pt=ak4jet_jerdown_p4[jetindexphoton12[1]]->Pt();
std::cout<<jet2_jerdown_pt<<std::endl;        
jet2_jerdown_eta=ak4jet_jerdown_p4[jetindexphoton12[1]]->Eta();
        jet2_jerdown_phi=ak4jet_jerdown_p4[jetindexphoton12[1]]->Phi();
        jet2_jerdown_e=ak4jet_jerdown_p4[jetindexphoton12[1]]->E();
               std::cout<<jet2_jerdown_pt<<" "<<jet2_jerdown_eta<<" "<<jet2_jerdown_phi<<" "<<jet2_jerdown_e<<std::endl;

        TLorentzVector j1_jerdown_p4;
        j1_jerdown_p4.SetPtEtaPhiE(jet1_jerdown_pt, jet1_jerdown_eta, jet1_jerdown_phi, jet1_jerdown_e);
        TLorentzVector j2_jerdown_p4;
        j2_jerdown_p4.SetPtEtaPhiE(jet2_jerdown_pt, jet2_jerdown_eta, jet2_jerdown_phi, jet2_jerdown_e);
        
        jet_jerdown_mass = (j1_jerdown_p4 + j2_jerdown_p4).M();
    }
}
    //***********************jer down jet****************************

     


     //*****************************************************************
   
       outTree_->Fill();
   }
   

//-------------------------------------------------------------------------------------------------------------------------------------//


void PKUTreeMaker::setDummyValues() {
     nVtx           = -1e1;
     triggerWeight  = -1e1;
     pileupWeight   = -1e1;
     lumiWeight     = -1e1;
     ptVlep         = -1e1;
     yVlep          = -1e1;
     phiVlep        = -1e1;
     massVlep       = -1e1;
     mtVlep         = -1e1;
     ptlep1         = -1e1;
     etalep1        = -1e1;
     philep1        = -1e1;
     met            = -1e1;
     metPhi         = -1e1;
     j1metPhi         = -1e1;
     j2metPhi         = -1e1;
    METraw_et = -99;
    METraw_phi = -99;
    METraw_sumEt = -99;
    MET_et = -99;
    MET_phi = -99;
    MET_sumEt = -99;
    MET_corrPx = -99;
    MET_corrPy = -99;
    
     lep            = -1e1;
 /*    ak4jet_pt[1] = -1e1;
     ak4jet_pt[2] = -1e1; 
     ak4jet_pt[3] = -1e1; 
     ak4jet_pt[4] = -1e1; 
     ak4jet_pt[5] = -1e1; 
     ak4jet_eta[0] = -1e1;
     ak4jet_eta[1] = -1e1;
     ak4jet_eta[2] = -1e1;
     ak4jet_eta[3] = -1e1;
     ak4jet_eta[4] = -1e1;
     ak4jet_eta[5] = -1e1;
     ak4jet_phi[0] = -1e1;
     ak4jet_phi[1] = -1e1;
     ak4jet_phi[2] = -1e1;
     ak4jet_phi[3] = -1e1;
     ak4jet_phi[4] = -1e1;
     ak4jet_phi[5] = -1e1;
     ak4jet_e[0] = -1e1;
     ak4jet_e[1] = -1e1;
     ak4jet_e[2] = -1e1;
     ak4jet_e[3] = -1e1;
     ak4jet_e[4] = -1e1;
     ak4jet_e[5] = -1e1;
     ak4jet_csv[0] = -1e3;
     ak4jet_csv[1] = -1e3;
     ak4jet_csv[2] = -1e3;
     ak4jet_csv[3] = -1e3;
     ak4jet_csv[4] = -1e3;
     ak4jet_csv[5] = -1e3;
     ak4jet_icsv[0] = -1e3;
     ak4jet_icsv[1] = -1e3;
     ak4jet_icsv[2] = -1e3;
     ak4jet_icsv[3] = -1e3;
     ak4jet_icsv[4] = -1e3;
     ak4jet_icsv[5] = -1e3;
  */

//     ak4jet_eta[6] = {-1e1, -1e1, -1e1, -1e1, -1e1, -1e1};
//     ak4jet_phi[6] = {-1e1, -1e1, -1e1, -1e1, -1e1, -1e1};
//     ak4jet_e[6] = {-1e1, -1e1, -1e1, -1e1, -1e1, -1e1};
//     ak4jet_csv[6] = {-1e1, -1e1, -1e1, -1e1, -1e1, -1e1};
//     ak4jet_icsv[6] = {-1e1, -1e1, -1e1, -1e1, -1e1, -1e1};

     photon_pt[0] = -1e1;
     photon_pt[1] = -1e1;
     photon_pt[2] = -1e1;
     photon_pt[3] = -1e1;
     photon_pt[4] = -1e1;
     photon_pt[5] = -1e1;
     photon_eta[0] = -1e1;
     photon_eta[1] = -1e1;
     photon_eta[2] = -1e1;
     photon_eta[3] = -1e1;
     photon_eta[4] = -1e1;
     photon_eta[5] = -1e1;
     photon_phi[0] = -1e1;
     photon_phi[1] = -1e1;
     photon_phi[2] = -1e1;
     photon_phi[3] = -1e1;
     photon_phi[4] = -1e1;
     photon_phi[5] = -1e1;
     photon_e[0] = -1e1;
     photon_e[1] = -1e1;
     photon_e[2] = -1e1;
     photon_e[3] = -1e1;
     photon_e[4] = -1e1;
     photon_e[5] = -1e1;

     photonet=-1e1;
     photoneta=-1e1;
     photonphi=-1e1;
     photone=-1e1;
     iphoton=-1;
     drla=1e1;

     jet1pt=-1e1;
     jet1eta=-1e1;
     jet1phi=-1e1;
     jet1e=-1e1;
     jet1csv=-1e1;
     jet1icsv=-1e1;
     jet2pt=-1e1;
     jet2eta=-1e1;
     jet2phi=-1e1;
     jet2e=-1e1;
     jet2csv=-1e1;
     jet2icsv=-1e1;
     drj1a=1e1;
     drj2a=1e1;
     drj1l=1e1;
     drj2l=1e1;
     Mjj=-1e1;
     deltaeta=-1e1;
     zepp=-1e1;
    
    jet1_jesup_pt=-1e1;
    jet1_jesup_e=-1e1;
    jet2_jesup_pt=-1e1;
    jet2_jesup_e=-1e1;
    jet_jesup_mass=-1e1;

    jet1_jesdown_pt=-1e1;
    jet1_jesdown_e=-1e1;
    jet2_jesdown_pt=-1e1;
    jet2_jesdown_e=-1e1;
    jet_jesdown_mass=-1e1;
    
    jet1_jerup_pt=-1e1;
    jet1_jerup_e=-1e1;
    jet2_jerup_pt=-1e1;
    jet2_jerup_e=-1e1;
    jet_jerup_mass=-1e1;

    jet1_jerdown_pt=-1e1;
    jet1_jerdown_e=-1e1;
    jet2_jerdown_pt=-1e1;
    jet2_jerdown_e=-1e1;
    jet_jerdown_mass=-1e1;

     ptel1          = -1e9;
     etaSC1         = -1e9;
     dEtaIn1        = -1e9;
     dPhiIn1        = -1e9;
     hOverE1        = -1e9;
     full5x5_sigma1 = -1e9;
     ooEmooP1       = -1e9;
     d01            = -1e9;
     dz1            = -1e9;
     relIso1        = -1e9;
     missingHits1   = -1e9; 
     passConVeto1   = -1e9;
     el1passID      = -1e9;
     el1tightID      = -1e9;
     el1mediumID      = -1e9;
     el1looseID      = -1e9;
     mu1tightID     = -1e9;
}

// ------------ method called once each job just before starting event loop  ------------
void 
PKUTreeMaker::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
PKUTreeMaker::endJob() {
  std::cout << "PKUTreeMaker endJob()..." << std::endl;
}

//define this as a plug-in
DEFINE_FWK_MODULE(PKUTreeMaker);
