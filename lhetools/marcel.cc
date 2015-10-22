// main71.cc is a part of the PYTHIA event generator.
// Copyright (C) 2013 Richard Corke.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

/*
 * Simple example of fastjet analysis. Roughly follows analysis of:
 * Throw di-jet events, 
 * reconstruct jets with FastJet,
 * calculate density in the core of jet
 */


#include "Pythia8/Pythia.h"
// This is the minimal interface needed to access FastJet. 
// A more sophisticated interface is demonstrated in main72.cc.
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
//#include "fastjet/contrib/ValenciaPlugin.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/Selector.hh"
#include "TNtuple.h"
#include "TFile.h"
#include "TLorentzVector.h"

using namespace Pythia8;
using namespace fastjet;
using namespace std;

/// set up a class to give standard (by default E-scheme)
 /// recombination, with additional tracking of flavour information in
 /// the user_index. 
 ///
 /// If you use this, you must explicitly set the user index to 0 for
 /// non-flavoured particles (the default value is -1);
 ///
 /// This will work for native algorithms, but not for all plugins
 typedef fastjet::JetDefinition::DefaultRecombiner DefRecomb;
 class FlavourRecombiner : public  DefRecomb {
 public:
   FlavourRecombiner(fastjet::RecombinationScheme recomb_scheme = fastjet::E_scheme) : 
     DefRecomb(recomb_scheme) {};
 
   virtual std::string description() const {return DefRecomb::description()
       +" (with user index addition)";}
 
   /// recombine pa and pb and put result into pab
   virtual void recombine(const fastjet::PseudoJet & pa, 
                          const fastjet::PseudoJet & pb, 
                          fastjet::PseudoJet & pab) const {
     DefRecomb::recombine(pa,pb,pab);
     pab.set_user_index(pa.user_index() + pb.user_index());
   }
 };

vector<PseudoJet> sorted_by_m(const vector<PseudoJet> & jets){
  vector<double> m(jets.size());
  for (size_t i = 0; i < jets.size(); i++) {m[i] = jets[i].m();}
  return objects_sorted_by_values(jets, m);
}

double modifJD(double pt1, double pt2, double DR){
	   return pt1*pt2*pow(DR, 4); 
	 }
double helicity_separation(PseudoJet j1, PseudoJet j2){
   PseudoJet jet1=j1;
  PseudoJet jet2=j2;
  PseudoJet rest;
  jet1.unboost(rest);
  jet2.unboost(rest);
  double p1x=jet1(0);
  double p2x=jet2(0);
  double p1y=jet1(1);
  double p2y=jet2(1);
  double p1z=jet1(2);
  double p2z=jet2(2);
  double p21=sqrt(p1x*p1x+p1y*p1y+p1z*p1z);
  double p22=sqrt(p2x*p2x+p2y*p2y+p2z*p2z);
  return (p1x*p2x+p1y*p2y+p1z*p2z)/(p21*p22);
}



 


int main() {
// Settings
int  nEvent = 100000; 
//for testing use less events
//int nEvent = 500;
bool doMPI  = true;

//count events which pass filters
 int event_0=0; //events which pass acceptance cuts
 int event_1=0; //events which pass top tagger
 int event_2=0; //events which pass top and H taggers
 int event_3=0; //events which also pass third b_tagger
 int event_all=0; // all events


  
// Generator
 Pythia pythia("/home/mello/physics_programs/Pythia_installation/xmldoc");


TNtuple *nt = new TNtuple("ntuple","ntuple","Tpt:Trap:Tphi:Tm:Te:Hpt:Hrap:Hphi:Hm:He");

 // TFile *f = new TFile("/home/mello/physics_programs/Pythia_LHE_FastJet/ttH.root","RECREATE");
  TFile *f = new TFile("/home/mello/physics_programs/Pythia_LHE_FastJet/ttbb.root","RECREATE");
 // TFile *f = new TFile("/home/mello/physics_programs/Pythia_LHE_FastJet/ttZ.root","RECREATE");
 // TFile *f = new TFile("/home/mello/physics_programs/Pythia_LHE_FastJet/ttjj.root","RECREATE");

pythia.readString("Beams:frameType = 4");
// read LHE input from MadGraph

//pythia.readString("Beams:LHEF =/home/mello/physics_programs/MG5_aMC_v2_2_1/NLO_ttH/Events/run_04/events.lhe ");
pythia.readString("Beams:LHEF =/home/mello/physics_programs/MG5_aMC_v2_2_1/BG_ttbb_NLO/Events/run_01/events.lhe ");
//pythia.readString("Beams:LHEF =/home/mello/physics_programs/MG5_aMC_v2_2_1/BG_ttZ_NLO/Events/run_01/events.lhe ");
//pythia.readString("Beams:LHEF =/home/mello/physics_programs/MG5_aMC_v2_2_1/BG_ttjj_LO/Events/run_02/events.lhe ");


// Initialisation, pp @ 14 TeV

pythia.init();


// Histograms
Hist dSigma1("dr distance", 100,0,0.4);



std::vector <fastjet::PseudoJet> fjInputs;
std::vector <fastjet::PseudoJet> fjPartonInputs;



bool firstEvent = true;



// Begin event loop. Generate event. Skip if error.
for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

if (!pythia.next()) continue;
 else{  cout<<"-------------------------------------------------"<<endl;}

 event_all++;
if (pythia.info.atEndOfFile()) break; 

   
    // Reset Fastjet input
    fjInputs.resize(0);
    fjPartonInputs.resize(0);
    
    // Keep track of missing ET
    Vec4 missingETvec;

    TLorentzVector t1;
    TLorentzVector t2;
    TLorentzVector H;
    TLorentzVector bZ;
    TLorentzVector bH;
    TLorentzVector pbH;
    vector<int> bjets;
    // Loop over event record to decide what to pass to FastJet
    
    for (int i = 0; i < pythia.event.size(); ++i) {
     

      //btag: put all bottoms' id in a vector 
      if(pythia.event[i].idAbs()==5){
	bjets.push_back(i);	
      }
     
     
     

      // Final state only

      // if ((!pythia.event[i].isFinal())) cout << "id = " << i << " " << pythia.event[i].idAbs() << " " << pythia.event[i].isFinal() << " " << pythia.event[i].daughter1() << " " << pythia.event[i].status() << endl;
      int didx = pythia.event[i].daughter1();
      
      
      if (!pythia.event[i].isFinal())        continue;
     
     



      //one hadronic and one leptonic top decay
      int midx = pythia.event[i].mother1();
      
      /*  //if first mother is a top
      if(pythia.event[midx].idAbs()==6 ||pythia.event[midx].idAbs()==-6){
	//get W index from decay
	int decay1=pythia.event[midx].daughter1();
	int decay2=pythia.event[midx].daughter2();

      }
      */ 


      // No neutrinos from Z for ZH -> vvbb analysis
      //One leptonic and one hadronic top decays
      /*   int midx = pythia.event[i].mother1();
      if (pythia.event[midx].idAbs() == 23 && 
	  (pythia.event[i].idAbs() == 12 || pythia.event[i].idAbs() == 14 ||
	   pythia.event[i].idAbs() == 16))     continue;
      */

      //no neutrinos releaved
      if (pythia.event[i].idAbs() == 12 || pythia.event[i].idAbs() == 14 ||
	   pythia.event[i].idAbs() == 16)     continue;

 
      // Store as input to Fastjet


      //check if the particle comes from a b (btag)
      int btag=0;
      for(int j=0; j< bjets.size(); j++){
	if(pythia.event[i].isAncestor(bjets[j])){btag=1; continue;cout<<"                       "<<j<<"       --------"<<endl;}
      }




      PseudoJet particle(pythia.event[i].px(), pythia.event[i].py(), pythia.event[i].pz(), pythia.event[i].e());
      particle.set_user_index(btag);
      
      fjInputs.push_back( particle );
      fjInputs[ fjInputs.size()-1].set_user_index(btag);
 
    }
    if (fjInputs.size() == 0 ) {
      cout << "Error: event with no final state particles" << endl;
      continue;
    }
 
  // Fastjet analysis - select algorithm and parameters
  double Rparam = 1.5;  //FATJET ANALYSIS
  FlavourRecombiner flav_recomb;    //it keeps track of  b-jets
  fastjet::JetAlgorithm           jetalgorithm=fastjet::cambridge_algorithm; 
  fastjet::Strategy               strategy = fastjet::Best;
  fastjet::RecombinationScheme    recombScheme = fastjet::E_scheme;
  fastjet::JetDefinition         *jetDef = NULL;
  //	if (Rparam > 1.55) Rparam=4;
  //jetDef = new fastjet::JetDefinition(fastjet::kt_algorithm,/*Rparam, */ recombScheme, strategy);
  
  jetDef = new fastjet::JetDefinition(jetalgorithm, Rparam, &flav_recomb,strategy);
 
  
  // Fastjet input
  
  
  // Run Fastjet algorithm
  
    vector <fastjet::PseudoJet> jt;
    //vector <fastjet::PseudoJet> jt2;
    vector <fastjet::PseudoJet> jH;
    vector <fastjet::PseudoJet> lep;
   
    
    fastjet::ClusterSequence clustSeq(fjInputs, *jetDef);
    
    // For the first event, print the FastJet details
    if (firstEvent) {
      cout << "Ran " << jetDef->description() << endl;
      cout << "Strategy adopted by FastJet was "
	   << clustSeq.strategy_string() << endl << endl;
      firstEvent = false;
    }
    
    vector <fastjet::PseudoJet>  fatj = sorted_by_pt(clustSeq.inclusive_jets());

    //   vector <fastjet::PseudoJet> jets = sorted_by_pt(clustSeq.inclusive_jets());
    //solo per far funzionare il cout.....è un check
    

  //check if there are at least two hard jet and a lepton
    
    /*  
     cout<<"size of the jet vector  "<<fatj.size()<<endl;
    
     cout <<   "        pt   rap   flav" << endl;
     for (unsigned i = 0; i < fatj.size(); i++) {
       cout << "jet " << i << ": "<< fatj[i].pt() << " "
	    << fatj[i].rap() << " " << fatj[i].user_index() << endl;}
    */
          
 
     
   


     vector<int> ntopcandidate;   //vector containg the jet numbers of top candidates
     vector<int> nHcandidate;    // vector containing the jet numbers of Higgs candidate
    //int extrabcandidate=0;  //extra b candidate
     vector<int> lepcandidate;   //leptoncandidate
     for(int i=0; i < fatj.size(); ++i){
       
       if(fatj[i].pt()>=200 && abs(fatj[i].rap())<=4 && fatj[i].rap()>2.5){
	 ntopcandidate.push_back(i);
	 
       }
       if(fatj[i].pt()>=200 && abs(fatj[i].rap())<=2.5 ){
	 ntopcandidate.push_back(i);
	 nHcandidate.push_back(i);
	 
	 
       }
       if(fatj[i].pt()>=15 && fatj[i].pt()<200 && fatj[i].rap()<=2.5){
	 lepcandidate.push_back(i);
       }
     }
     if(ntopcandidate.size()==0 || nHcandidate.size()==0 ||  lepcandidate.size()==0 ||(nHcandidate.size()==1 && ntopcandidate.size()==1) ){
       cout<<"Not enough fat hard jets!!"<<endl<<"Two hard jets and a lepton are required "<<endl;
       event_0++;
       continue;
      
     }
     // cout<<"number of fat jets candidates "<<ntopcandidate.size()<<"  (of whom Higgs= "<<nHcandidate.size()<<" )"<<endl; 

 
    
     //TOP-TAGGER
     
     //----------------------------------------------------------------------
     
     //counts the number of substructures for each fat jet i
     //initialization
     vector<int> aux;
     for(int l=0; l<ntopcandidate.size();l++){
       aux.push_back(0);
     }
     
     
     vector<PseudoJet> constituents;
     int nsubstructure=0;  //total number of substructures found
     for(int l=0; l<ntopcandidate.size(); l++){
       
       
       constituents.push_back(fatj[ntopcandidate[l]]);
       //cout<<"pt and mass of fatjet passed to subroutine "<< constituents[constituents.size()-1].pt()  <<"    "<<  constituents[constituents.size()-1].m() <<endl;
      
       nsubstructure=constituents.size()-1;
       
       
       // issue a warning if the jet is not obtained through a C/A
       // clustering
       /*     if ((! fatj[ntopcandidate[i]].has_associated_cluster_sequence()) ||
	      (fatj[ntopcandidate[i]].validated_cs()->jet_def().jet_algorithm() != cambridge_algorithm))
	      cout<<"MassDropTagger should only be applied on jets from a Cambridge/Aachen clustering; use it with other algorithms at your own risk."<<endl;
       */
       
       PseudoJet j1, j2;
       bool had_parents; 
      
     
       while ((had_parents=constituents[nsubstructure].has_parents(j1, j2)) &&  constituents[nsubstructure].m()>=30) {
	 int divide=0;
	
	 //cout<< "        nsubstructure "<<nsubstructure<<"                       constituents size "<<constituents.size()<<endl;

	 // make parent1 the more massive jet
	 if (j1.m2() < j2.m2()) swap(j1,j2);
	 // if we pass the conditions on the mass drop, then add the more
	 // massive to the list of interesting substructers, otherwise add both.
	 
	 if ( (j1.m() >  constituents[nsubstructure].m()*0.8)) {
	   constituents[nsubstructure]=j1;
	  
	   
	 }
	 else{
	   
	   constituents[nsubstructure]=j1;
     
	   constituents.push_back(j2);
	   divide=1;
	 }
	 if (constituents[nsubstructure].m()<30 && divide==1){
	   nsubstructure++;
	  
	   //if(constituents[nsubstructure].m()<30) continue; //not needed
	 }
	 if (!had_parents) {cout<<"Stable particle with m>30Gev!"<<endl; continue;}//is it possible or error?
	
       }
       
       int sum=0;
       if(l==0)   aux[l]=constituents.size();
       else{
	 for(int j=0;j<l;j++) sum=sum+aux[j];
	 aux[l]=constituents.size()-sum;
	
       }
     
       /*
       cout<<"number of relevant substructures of "<<l<<" : "<<aux[l]<<endl;
      
       for (unsigned j = 0; j < aux[l]; j++) {
	 cout << "    constituent " << j << "’s mass: "<< constituents[constituents.size()-aux[l]+j].m() << endl;
       }
       */
     
  


       //here put filter

        
       double   Rfilt =0.3; // arbitrary... (look for references)
       unsigned nfilt = 3;               // number of pieces we'll take
       
       Filter filter(JetDefinition(cambridge_algorithm, Rfilt, &flav_recomb),
		     SelectorNHardest(nfilt));
     
       for(int i=sum;i<sum+aux[l];i++){
	 constituents[i]=filter(constituents[i]);
       }
       

       vector<PseudoJet> Wbosons;
       vector<PseudoJet> Top;
       vector<int> associatedW;
       for(int i =sum; i<sum+aux[l]; i++){
	 for(int j=i+1; j<sum+aux[l];j++){
	   PseudoJet j3;
	   jetDef->recombiner()->recombine(constituents[j], constituents[i], j3);
	  

	   if(j3.m()>=65 && j3.m()<=95) {
	     Wbosons.push_back(j3);



	     //here put second filter
	      
	       for(int k=sum;k<sum+aux[l];k++){
	       if(k!=i && k!=j){
		 constituents[k]=filter(constituents[k]);
	       }
	       if(k==i || k== j)	 constituents[k]=filter(constituents[k]);;
	     }
	     
	     for(int k=sum;k<sum+aux[l];k++){
	       if(k!=i && k!=j){
		 PseudoJet j4=j3+constituents[k];
		 
		   
		 if (j4.m()>=150 && j4.m()<=200 && helicity_separation(j3, j4)<0.7){
		   cout<<"helicity separation: "<<helicity_separation(j3, j4)<<endl;		   	
		  
		   associatedW.push_back(Wbosons.size()-1);
		   Top.push_back(j4);
		 }
	       }
	     }
	   }
	 }
       }
       
       if(Wbosons.size()==0 ||Top.size()==0){
	 PseudoJet Tout;
	 jt.push_back(Tout); //add an empty jet if there is no top
       }
       else{
	 //if more than one top, choose the one that minimize choosebestmass
	 PseudoJet Tout;
	 double choosebestmass=10000000000;
	 for(int i=0;i<Top.size();i++){       
	   double var=abs(Top[i].m()-173.2)+abs(Wbosons[associatedW[i]].m()-80.2);
	   if(var<choosebestmass){
	     choosebestmass=var;
	     Tout=Top[i];
	   }
	 }  
	 jt.push_back(Tout); //result of top tagger
       }
     }

     
     if(jt.size()==0){cout<< "NO TOP JETS"<<endl; event_1++; continue;}
     

   
     
     for(int i=0; i<jt.size();i++){
       // cout<<"Top output "<<i<<" with mass "<<jt[i].m()<<endl;
     }
     
     //if more than a top, choose the one with mass nearer to top pole mass
     PseudoJet bestT;
     int fatjetid;
     if(jt.size()>=1){
       
       double choosebestmass=173.2; // a mass 0 candidate is rejected
       
       for(int i=0;i<jt.size();i++){       
	 double var=abs(jt[i].m()-173.2);
	 if(var<=choosebestmass){
	   choosebestmass=var;
	   bestT=jt[i];
	   fatjetid=i;
	 }
       }
       if(fatjetid>jt.size()) {cout<<"ERROR1"<<endl;return 0; }
     }// else{fatjetid=0;}
     
     if(bestT.m()==0){cout<<"NO TOP CANDIDATES!"<<endl; event_1++; continue;}
     else{
     cout<<endl<<"Best top is from fatjet "<<fatjetid<<" with mass "<<bestT.m()<<endl;;
     }
    



     //H-TAGGER
     
     //----------------------------------
    
          //counts the number of substructures for each fat jet i
     //initialization
     vector<int> auxH;
     
     for(int l=0; l<nHcandidate.size();l++){
       auxH.push_back(0);
     }
     
     
     vector<PseudoJet> constituentsH;
     int nsubstructureH=0;  //total number of substructures found
     for(int l=0; l<nHcandidate.size(); l++){
      
       if(nHcandidate[l]!=fatjetid){
	 
	 constituentsH.push_back(fatj[nHcandidate[l]]);
	 constituentsH[constituentsH.size()-1].set_user_index(fatj[nHcandidate[l]].user_index());
	 // cout<<"pt and flavour of fatjet passed to subroutine "<< constituentsH[constituentsH.size()-1].pt()  <<"    "<<  constituentsH[constituentsH.size()-1].user_index() <<endl;
	 
	 nsubstructureH=constituentsH.size()-1;
	 
	 
	 // issue a warning if the jet is not obtained through a C/A
	 // clustering
	 /*     if ((! fatj[ntopcandidate[i]].has_associated_cluster_sequence()) ||
		(fatj[ntopcandidate[i]].validated_cs()->jet_def().jet_algorithm() != cambridge_algorithm))
		cout<<"MassDropTagger should only be applied on jets from a Cambridge/Aachen clustering; use it with other algorithms at your own risk."<<endl;
	 */
	 
	 PseudoJet j1, j2;
	 bool had_parents; 
	  
	 
	 while ((had_parents=constituentsH[nsubstructureH].has_parents(j1, j2)) &&  constituentsH[nsubstructureH].m()>40) {
	   
	   
	   int divide=0;
	   // make parent1 the more massive jet
	   if (j1.m2() < j2.m2()) swap(j1,j2);
	   // if we pass the conditions on the mass drop, then add the more
	   // massive to the list of interesting substructers, otherwise add both.
	   
	   if ( (j1.m() >  constituentsH[nsubstructureH].m()*0.9)) {
	     constituentsH[nsubstructureH]=j1;
	     constituentsH[nsubstructureH].set_user_index(j1.user_index());
	     
	   }
	   else{
	     constituentsH[nsubstructureH]=j1;
	     constituentsH[nsubstructureH].set_user_index(j1.user_index());
	     constituentsH.push_back(j2);
	     divide=1;
	   }
	   if (constituentsH[nsubstructureH].m()<40 && divide==1){
	     nsubstructureH++;
	     if(constituentsH[nsubstructureH].m()<30) continue;
	   }
	   if (!had_parents) cout<<"Stable particle with m>40Gev!"<<endl;//is it possible or error?	
	 }
	
	 int sum=0;
	 if(l==0)   auxH[l]=constituentsH.size();
       else{
	 for(int j=0;j<l;j++) sum=sum+auxH[j];
	 auxH[l]=constituentsH.size()-sum;
	 
       }
	 
	 /*
	   cout<<"number of relevant substructures of "<<l<<" : "<<aux[l]<<endl;
	   
	   for (unsigned j = 0; j < aux[l]; j++) {
	   cout << "    constituent " << j << "’s mass: "<< constituents[constituents.size()-aux[l]+j].m() << endl;
	   }
	 */
	 
	 
	 //filter
	 
	 double   Rfilt =0.3; // arbitrary... (look for references)
	 unsigned nfilt = 3;               // number of pieces we'll take
	 
	 Filter filter(JetDefinition(cambridge_algorithm, Rfilt, &flav_recomb),
		       SelectorNHardest(nfilt));
	 
	 
	 for(int i=sum;i<sum+auxH[l];i++){
	   constituentsH[i]=filter(constituentsH[i]);
	 }
	 

	 vector<PseudoJet> B;
	 vector<PseudoJet> Higgs;
	 // vector<PseudoJet> pairs;
	 //vector<int> b1id;
     
	  vector<double> dist;
	
	 for(int i =sum; i<sum+auxH[l]; i++){
	   for(int j=i+1; j<sum+auxH[l];j++){
	     PseudoJet j3;
	     //se due jet sono da bottom e antibottom ricostruisce higgs
	     if(constituentsH[i].user_index()!=0 && constituentsH[j].user_index()!=0){
	       j3+=constituentsH[i];
	       j3+=constituentsH[j];
	      
	       if(j3.m()>115 && j3.m()<135){Higgs.push_back(j3);} //mass window
	       //Higgs.push_back(j3);   //without mass window
	     }
	     /* double Jdist= modifJD(constituentsH[i].perp(),constituentsH[i].perp(),constituentsH[i].delta_R(constituentsH[j]));
	     
	     j3+=constituentsH[i];
	     j3+=constituentsH[j];
	     pairs.push_back(j3);
	     dist.push_back(Jdist);
	     */
	   }
	 }
	 // objects_sorted_by_values(pairs, dist);//ordered pairs by modified Jade distance
	 
	 //here filter the first three leading pairs
	    
	   
	 
       
    

	 
	 if(Higgs.size()==0){
	   PseudoJet Hout;
	   jH.push_back(Hout); //add an empty jet if there is no higgs
	 }
	 else{
	   //if more than one higgs, choose the one that minimize choosebestmass
	   PseudoJet Hout;
	   double choosebestmass=10000000000;
	   for(int i=0;i<Higgs.size();i++){       
	     double var=abs(Higgs[i].m()-125);
	     if(var<=choosebestmass){
	       choosebestmass=var;
	       Hout=Higgs[i];
	     }
	   }
	   jH.push_back(Hout); //result of higgs tagger
	   
	 }
	 
       }
       else{PseudoJet Hout; jH.push_back(Hout);} //for the jet already used generate a M=0 Higgs (it will be discarded) 
      
     }


     if(jH.size()==0){cout<<"NO HIGGS CANDIDATE!"<<endl;event_2++;continue;}
     for(int i=0; i<jH.size();i++){
       //  cout<<"Higgs output "<<i<<" with mass "<<jH[i].m()<<endl;
     }
     
     PseudoJet bestH;
     int fatjetidH=0;

     if(jH.size()>=1){
       
       double choosebestmass=125.2; // a mass0 candidate is rejected
       
       for(int i=0;i<jH.size();i++){       
	 double var=abs(jH[i].m()-125.2);
	 if(var<choosebestmass){
	   choosebestmass=var;
	   bestH=jH[i];
	   fatjetidH=i;
	 }
       }
      
       // if(fatjetidH>jH.size()) {cout<<"ERROR"<<endl; return 0;}
     }
     if(bestH.m()==0) {cout<<"NO HIGGS CANDIDATES!"<<endl;event_2++; continue;}
     cout<<endl<<"Best Higgs is from fatjet "<<fatjetidH<<" with mass "<<bestH.m()<<endl;
	


     ///THIRD B-TAG
     //---------------------------------------------------------------

     int third_btag=0;
     vector<PseudoJet> remaining;
     for(int i=0; i < fatj.size(); ++i){
       if(i!=fatjetid && i!=fatjetidH) remaining.push_back(fatj[i]);
     }

     //recluster remaining particles with R=0.6
      fastjet::JetDefinition         *jetDef2 = NULL;
 
      jetDef2 = new fastjet::JetDefinition(jetalgorithm, 0.6, &flav_recomb,strategy);
      ClusterSequence reclust(remaining, *jetDef2);
      for(int i=0; i<remaining.size(); i++){
	if(remaining[i].perp()>30 && 
	   abs(remaining[i].rap())<2.5 &&
	   remaining[i].delta_R(fatj[fatjetid])>0.4 &&
	   remaining[i].delta_R(fatj[fatjetidH])>0.4 &&
	   remaining[i].user_index()!=0 ){
	  third_btag=1;
	 
	}
      }
      if(third_btag==0)event_3++;
 
      //count events after every cut, jet tag, ecc...
      
     
      /* if(bestT.m()!=0 && jt.size()!=0) event_1++;
      if(bestT.m()!=0 && bestH.m()!=0 &&  jt.size()!=0 && jH.size()!=0 ) event_2++;
      */
     //if there are H and t and the third_btag, fill histo
     if(bestT.m()!=0 && bestH.m()!=0 &&  jt.size()!=0 && jH.size()!=0 && third_btag!=0){
       //event_3++;
       cout<<endl;
       t1.SetPxPyPzE(
						  bestT.px(),
						  bestT.py(), 
						  bestT.pz(), 
						  bestT.e());
       H.SetPxPyPzE(
						  bestH.px(),
						  bestH.py(), 
						  bestH.pz(), 
						  bestH.e());
       vector<float> input;
       input.push_back(bestT.perp());
       input.push_back(bestT.rap());
       input.push_back(bestT.phi_std());
       input.push_back(bestT.m());
       input.push_back(bestT.e());

       input.push_back(bestH.perp());
       input.push_back(bestH.rap());
       input.push_back(bestH.phi_std());
       input.push_back(bestH.m());
       input.push_back(bestH.e());
       nt->Fill(&input[0]);
     }



 } 
 cout<<"All events:  "<<event_all<<endl;
 cout<<"Events after acceptance cuts: "<<event_all-event_0<<endl;
 cout<<"Events with one top tag: "<<event_all-event_0-event_1<<endl;
 cout<<"Events with Higgs tag: "<<event_all-event_0-event_1-event_2<<endl;
 cout<<"Events with third b tag: "<<event_all-event_0-event_1-event_2-event_3<<endl;
 nt->Write();
 f->Write();
 f->Close();
     // Statistics
     //pythia.stat();
     /*  
     
     ntav->Write();
     nt->Write();
     f->Write();
     f->Close();
     */
     
     
     // Done.
     
     return 0;
}

