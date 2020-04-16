void ttbar750(){
	#include <iostream>
	#include <string>
	#include <stdio.h> 


	//SIGNAL PART

	// Returns the new jet flavour, 5 = b, 4 = c, 15 = tau, 0 = light
	TFile *file = TFile::Open("/Users/sinaaktas/DERSLER/PHYS492/ttData/mc_110902.ZPrime750.root");
	TTree *tree = (TTree*) file->Get("mini");

	const int vs =20;

	UInt_t jet_n; //number of selected jets

	Int_t //vector<int>
			jet_trueflav[vs];  //true flavor of the jet

	double met_et; //Transverse energy of the missing momentum vector

	Float_t //vector<float>
       		jet_pt[vs],   //transverse momentum 
    	    jet_eta[vs],  //pseudorapidity
       		jet_phi[vs],  //azimuthal angle
        	jet_E[vs],    //energy
        	jet_m[vs];    //invariant mass of the jet


	tree->SetBranchAddress("jet_n", &jet_n);
	tree->SetBranchAddress("jet_trueflav", &jet_trueflav);
	tree->SetBranchAddress("met_et", &met_et);
	tree->SetBranchAddress("jet_pt", &jet_pt);
	tree->SetBranchAddress("jet_eta", &jet_eta);
	tree->SetBranchAddress("jet_phi", &jet_phi);
	tree->SetBranchAddress("jet_E", &jet_E);
	tree->SetBranchAddress("jet_m", &jet_m);


	//Arrange the histograms
	TH1F *W1hist = new TH1F("W1hist","Hadronic W Boson recovery(GeV)",40,0,200);
	W1hist->SetFillColor(kBlue);

	TH1F *W2hist = new TH1F("W2hist","Hadronic W Boson recovery(GeV)",40,0,200);
	W2hist->SetFillColor(kBlue);

	TH1F *t1hist = new TH1F("t1hist","Hadronic Top quark recovery(GeV)",40,0,700);
	t1hist->SetFillColor(kRed);

	TH1F *t2hist = new TH1F("t2hist","Hadronic Top quark recovery(GeV)",40,0,700);
	t2hist->SetFillColor(kRed);

	TH1F *Zhist = new TH1F("Zhist","Hadronic Z' recovery(GeV)",40,0,2000);
	Zhist->SetFillColor(kGreen);


	// IMPORTANT: fraction events we want to run
	int nentries;
	nentries = (Int_t)tree->GetEntries();
	Float_t fraction_events = 1.; // fraction of all events
	Float_t events_to_run = nentries*fraction_events; // number of events we run


	std::cout << "Total # events = "  << nentries
	          << ". Events to run = " << events_to_run
	          << " corresponding to " << fraction_events*100
	          << "% of total events!" << std::endl;   // output of number of events run 

	printf("hello");



	Bool_t  bool_E_240; // The trigger for the Top Template Tagger selection requires an event to have at least one anti-kt jet with a distance parameter R = 1.0 and ET > 240 GeV.
	
	Int_t 
		lightjet[10] , // list of order of light jet in the one event jet list
		bjet[10];      // list of order of   b   jet in the one event jet list

	int 
		nbytes, //to select event
		bjet_n, //number of b jet in one event
		lightjet_n; // number of light jet in one event



	Float_t 
			wchi2,   //chi2 for W mass
			topchi2, //chi2 for top mass
			chi2check = 9999.;

    // Define Lorentz Vector of jets, W+ , W- , t ,tbar, Z'
	TLorentzVector j1;
	TLorentzVector j2;
	TLorentzVector j3;
	TLorentzVector j4;
	TLorentzVector b1;
	TLorentzVector b2;

	TLorentzVector W1invmass;
	TLorentzVector W2invmass;
	TLorentzVector T1invmass;
	TLorentzVector T2invmass;
	TLorentzVector Zinvmass;






	for (Int_t i=0; i<events_to_run; i++){ // begin for loop in for each event
		nbytes = tree->GetEntry(i);        // chose event number
		if(jet_n >5 && met_et < 200000.){  // to produce Z' we need at least 6 jet, (it's a pre-selection)
			bool_E_240=0; //set 0 to boolean of energy of at least one jet >240 GeV 
			bjet_n=0;     // set 0 to number of   b   jet in the event
			lightjet_n=0; // set 0 to number of light jet in the event

			for(Int_t t=0; t<15; t++){  // set zero to all variable in both b jet and light jet list
				lightjet[t]=0;
				bjet[t]=0;

			}

			for(Int_t h=0; h<jet_n; h++){  // for loop for all jet in one event
				if(jet_E[h]>240000.){   // check that if any one has energy more than 240 GeV
					bool_E_240 = 1;     //make boolean True
				}
				if(jet_pt[h]> 25000. && TMath::Abs(jet_eta[h]) < 2.5 && jet_trueflav[h]==5)	{  //we select b jet has transverse momentum more than 25 GeV, jet_trueflav==5 is mean that it is b jet
					bjet[bjet_n]=h; // set position of b jet into b jet list
					bjet_n++;       // increase number of b jet +1
				}
				else if(jet_pt[h]> 25000. && jet_trueflav[h]==0){ //we select light jet has transverse momentum more than 25 GeV, jet_trueflav==0 is mean that it is light jet
					lightjet[lightjet_n]=h;  // set position of b jet into light jet list
					lightjet_n++;    		 // increase number of light jet +1
				}
			}
			if(bool_E_240 == 1 && bjet_n >= 2 && lightjet_n >=4){ // we need at least 2 b jet and 4 light jet

				// In this part we pair all combination of light jet to produce W+ and W-
				// to decide goodJet; we check their chi2 (True value of W Boson Mass is 80.385±0.015 GeV/c2)
				for(Int_t j_1=0; j_1<lightjet_n; j_1++){
					for(Int_t j_2=0; j_2<lightjet_n; j_2++){
						for(Int_t j_3=0; j_3<lightjet_n; j_3++){
							for(Int_t j_4=0; j_4<lightjet_n; j_4++){
								if(j_1!=j_2 && j_1!=j_3 && j_1!=j_4 && j_2!=j_3 && j_2!=j_4 && j_3!=j_4){ //to prevent duplicate jets

									j1.SetPtEtaPhiE(jet_pt[lightjet[j_1]], jet_eta[lightjet[j_1]], jet_phi[lightjet[j_1]], jet_E[lightjet[j_1]]);
           							j2.SetPtEtaPhiE(jet_pt[lightjet[j_2]], jet_eta[lightjet[j_2]], jet_phi[lightjet[j_2]], jet_E[lightjet[j_2]]);
           							j3.SetPtEtaPhiE(jet_pt[lightjet[j_3]], jet_eta[lightjet[j_3]], jet_phi[lightjet[j_3]], jet_E[lightjet[j_3]]);
           							j4.SetPtEtaPhiE(jet_pt[lightjet[j_4]], jet_eta[lightjet[j_4]], jet_phi[lightjet[j_4]], jet_E[lightjet[j_4]]);

           							W1invmass = j1 +j2;
									W2invmass = j3 +j4;

									// In this part we pair all combination of b jet to produce t and tbar with W bosons we found previous part
									// to decide goodJet; we check their chi2
									for(Int_t b_1=0; b_1<lightjet_n; b_1++){
										for(Int_t b_2=0; b_2<bjet_n; b_2++){
											if(b_1!=b_2){ //to prevent duplicate jets


											
												b1.SetPtEtaPhiE(jet_pt[bjet[b_1]], jet_eta[bjet[b_1]], jet_phi[bjet[b_1]], jet_E[bjet[b_1]]);
           										b2.SetPtEtaPhiE(jet_pt[bjet[b_2]], jet_eta[bjet[b_2]], jet_phi[bjet[b_2]], jet_E[bjet[b_2]]);

           										T1invmass = W1invmass + b1;
												T2invmass = W2invmass + b2;

												Float_t topchi2;
												Float_t wchi2; 
												topchi2= TMath::Power(((T1invmass.M()/1000.)-(T2invmass.M()/1000.))/4.2 , 2.);
												wchi2 = TMath::Power(((W1invmass.M()/1000.)-80.4)/2.1, 2.) + TMath::Power(((W2invmass.M()/1000.)-80.4)/2.1 , 2.);


												// we check the chi2 is different from other chi2 of jet pair in same event to prevent count same Z' more than one times
												// we check small chi2 value
												// we check transverse momentum of top quarks more than 175 GeV
												if ( TMath::Abs(topchi2+wchi2) !=chi2check && TMath::Abs(topchi2+wchi2)<80. && T1invmass.Pt()>200000. && T2invmass.Pt()>200000. ){ 

													chi2check = TMath::Abs(topchi2+wchi2); // set chi2 value for compare next jet pair

													Zinvmass = T1invmass + T2invmass;  //combine t and t bar lorentz vector to find Z'

													//Fill the histograms
													W1hist->Fill(   (W1invmass.M())/1000.      );
													W2hist->Fill(   (W2invmass.M())/1000.      );

													t1hist->Fill(   (T1invmass.M())/1000.      );
													t2hist->Fill(   (T2invmass.M())/1000.      );

													Zhist->Fill(   (Zinvmass.M())/1000.      );
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}



	//BACKGROUND PART (same code as signal part without ttbar quark have transverse momentum higher than 200 GeV)

	TFile *BGfile = TFile::Open("/Users/sinaaktas/DERSLER/PHYS492/ttData/mc_117049.ttbar_had.root");
	TTree *BGtree = (TTree*) BGfile->Get("mini");


	TH1F *Backgroundhist = new TH1F("Backgroundhist","Z' with Background (GeV)",40,0,2000);
	Backgroundhist->SetFillColor(kMagenta);


	// IMPORTANT: fraction events we want to run
	int 
		BGnbytes,
		BGnentries;

	BGnentries = (Int_t)BGtree->GetEntries();

	Float_t fraction_BGevents = 1.; // fraction of all events
	Float_t BGevents_to_run = BGnentries*fraction_BGevents; // number of events we run

	BGtree->SetBranchAddress("jet_n", &jet_n);
	BGtree->SetBranchAddress("jet_trueflav", &jet_trueflav);
	BGtree->SetBranchAddress("met_et", &met_et);
	BGtree->SetBranchAddress("jet_pt", &jet_pt);
	BGtree->SetBranchAddress("jet_eta", &jet_eta);
	BGtree->SetBranchAddress("jet_phi", &jet_phi);
	BGtree->SetBranchAddress("jet_E", &jet_E);
	BGtree->SetBranchAddress("jet_m", &jet_m);


	for (Int_t i=0; i<BGevents_to_run; i++){ // begin for loop in for each event
		BGnbytes = BGtree->GetEntry(i);        // chose event number
		if(jet_n >5 && met_et < 200000.){  // to produce Z' we need at least 6 jet, (it's a pre-selection)
			bool_E_240=0; //set 0 to boolean of energy of at least one jet >240 GeV 
			bjet_n=0;     // set 0 to number of   b   jet in the event
			lightjet_n=0; // set 0 to number of light jet in the event

			for(Int_t t=0; t<15; t++){  // set zero to all variable in both b jet and light jet list
				lightjet[t]=0;
				bjet[t]=0;

			}

			for(Int_t h=0; h<jet_n; h++){  // for loop for all jet in one event
				if(jet_E[h]>240000.){   // check that if any one has energy more than 240 GeV
					bool_E_240 = 1;     //make boolean True
				}
				if(jet_pt[h]> 25000. && TMath::Abs(jet_eta[h]) < 2.5 && jet_trueflav[h]==5)	{  //we select b jet has transverse momentum more than 25 GeV, jet_trueflav==5 is mean that it is b jet
					bjet[bjet_n]=h; // set position of b jet into b jet list
					bjet_n++;       // increase number of b jet +1
				}
				else if(jet_pt[h]> 25000. && jet_trueflav[h]==0){ //we select light jet has transverse momentum more than 25 GeV, jet_trueflav==0 is mean that it is light jet
					lightjet[lightjet_n]=h;  // set position of b jet into light jet list
					lightjet_n++;    		 // increase number of light jet +1
				}
			}
			if(bool_E_240 == 1 && bjet_n >= 2 && lightjet_n >=4){ // we need at least 2 b jet and 4 light jet

				// In this part we pair all combination of light jet to produce W+ and W-
				// to decide goodJet; we check their chi2 (True value of W Boson Mass is 80.385±0.015 GeV/c2)
				for(Int_t j_1=0; j_1<lightjet_n; j_1++){
					for(Int_t j_2=0; j_2<lightjet_n; j_2++){
						for(Int_t j_3=0; j_3<lightjet_n; j_3++){
							for(Int_t j_4=0; j_4<lightjet_n; j_4++){
								if(j_1!=j_2 && j_1!=j_3 && j_1!=j_4 && j_2!=j_3 && j_2!=j_4 && j_3!=j_4){ //to prevent duplicate jets

									j1.SetPtEtaPhiE(jet_pt[lightjet[j_1]], jet_eta[lightjet[j_1]], jet_phi[lightjet[j_1]], jet_E[lightjet[j_1]]);
           							j2.SetPtEtaPhiE(jet_pt[lightjet[j_2]], jet_eta[lightjet[j_2]], jet_phi[lightjet[j_2]], jet_E[lightjet[j_2]]);
           							j3.SetPtEtaPhiE(jet_pt[lightjet[j_3]], jet_eta[lightjet[j_3]], jet_phi[lightjet[j_3]], jet_E[lightjet[j_3]]);
           							j4.SetPtEtaPhiE(jet_pt[lightjet[j_4]], jet_eta[lightjet[j_4]], jet_phi[lightjet[j_4]], jet_E[lightjet[j_4]]);

           							W1invmass = j1 +j2;
									W2invmass = j3 +j4;

									// In this part we pair all combination of b jet to produce t and tbar with W bosons we found previous part
									// to decide goodJet; we check their chi2
									for(Int_t b_1=0; b_1<lightjet_n; b_1++){
										for(Int_t b_2=0; b_2<bjet_n; b_2++){
											if(b_1!=b_2){ //to prevent duplicate jets


											
												b1.SetPtEtaPhiE(jet_pt[bjet[b_1]], jet_eta[bjet[b_1]], jet_phi[bjet[b_1]], jet_E[bjet[b_1]]);
           										b2.SetPtEtaPhiE(jet_pt[bjet[b_2]], jet_eta[bjet[b_2]], jet_phi[bjet[b_2]], jet_E[bjet[b_2]]);

           										T1invmass = W1invmass + b1;
												T2invmass = W2invmass + b2;

												Float_t topchi2;
												Float_t wchi2; 
												topchi2= TMath::Power(((T1invmass.M()/1000.)-(T2invmass.M()/1000.))/4.2 , 2.);
												wchi2 = TMath::Power(((W1invmass.M()/1000.)-80.4)/2.1, 2.) + TMath::Power(((W2invmass.M()/1000.)-80.4)/2.1 , 2.);


												// we check the chi2 is different from other chi2 of jet pair in same event to prevent count same Z' more than one times
												// we check small chi2 value
												// we check transverse momentum of top quarks more than 175 GeV
												if ( TMath::Abs(topchi2+wchi2) !=chi2check && TMath::Abs(topchi2+wchi2)<80.){ 

													chi2check = TMath::Abs(topchi2+wchi2); // set chi2 value for compare next jet pair

													Zinvmass = T1invmass + T2invmass;  //combine t and t bar lorentz vector

													//Fill the histogram

													Backgroundhist->Fill(   (Zinvmass.M())/1000.      );
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}






	

	gStyle->SetTextFont(42);
   	auto cz = new TCanvas("cz", "", 800, 700);
   	cz->Divide(2,3);

	cz->cd(1);
	W1hist->Draw();
	cz->cd(2);
	W2hist->Draw();

	cz->cd(3);
	t1hist->Draw();
	cz->cd(4);
	t2hist->Draw();

	cz->cd(5);
	Zhist->Draw();
	cz->cd(6);
	Backgroundhist->Draw();
	Zhist->Draw("SAME");

	auto legend = new TLegend(0.7,0.6,1.,0.7);
	legend->AddEntry(Backgroundhist,"t+tbar-->Jets","f");
	legend->AddEntry(Zhist,"Z'  (MC)","f");
	legend->Draw();

	

	cz->Draw();







}