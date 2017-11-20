#include "TROOTAnalysis.hh"

TROOTAnalysis::TROOTAnalysis(){
}


TROOTAnalysis::TROOTAnalysis(std::unique_ptr<TChain> &ch) :

        EcalTree(ch->GetTree()),

        D(new TCanvas("DC","Direction vector difference")),
        dx1(new TH1D("DeltaX_photon1","DeltaX_photon1", 1000,-3,3)),
        dy1(new TH1D("DeltaY_photon1","DeltaY_photon1", 1000,-3,3)),
        dz1(new TH1D("DeltaZ_photon1","DeltaZ_photon1", 1000,-3,3)),
        dx2(new TH1D("DeltaX_photon1_noPCA","DeltaX_photon1_noPCA", 1000,-3,3)),
        dy2(new TH1D("DeltaY_photon1_noPCA","DeltaY_photon1_noPCA", 1000,-3,3)),
        dz2(new TH1D("DeltaZ_photon1_noPCA","DeltaZ_photon1_noPCA", 1000,-3,3)),
        projectionC(new TCanvas("ProjectionC", "ProjectionC")),
        projection_pca(new TH2D("Prjection", "Projection", 10000,-2000,0,10000,-2000,0)),
        projection_minimization(new TH2D("Prjection", "Projection", 10000,-2000,0,10000,-2000,0))

{

        std::cout << "assigned Tree" << '\n';

        nofEntries=EcalTree->GetEntries();

        Cevent=new B4ROOTEvent();
        EcalTree->SetBranchAddress("EventBranch", &Cevent);

        EcalTree->GetEntry(0);

        GapThickness=Cevent->GapThickness();
        AbsoThickness=Cevent->AbsoThickness();
        tiledimX=Cevent->TilesizeX();
        tiledimY=Cevent->TilesizeY();
        calsizeXY=Cevent->calsizeXY();
        nofLayers=Cevent->NumberOfLayers();
        absofirst=Cevent->GetAbsFirst();

        std::cout<<"Absorber thickness is "<<AbsoThickness<<"mm"<<std::endl;
        std::cout<<"Scintillator thickness is "<<GapThickness<<"mm"<<std::endl;
        std::cout<<"Tilesize is "<<tiledimX<<"*"<<tiledimY<<"mm^2"<<std::endl;
        std::cout << std::boolalpha;
        std::cout<<"AbsorberFirst: "<<absofirst<<std::endl;

        histsizeX=calsizeXY/tiledimX;
        histsizeY=calsizeXY/tiledimY;
        histsizeZ=nofLayers;

        std::unique_ptr<TH2D> _h1(new TH2D("h1", "h1", histsizeX+1,-0.5,histsizeX+0.5,histsizeY+1,-0.5,histsizeY+0.5));
        std::unique_ptr<TH2D> _h2(new TH2D("h2", "h2", histsizeX+1,-0.5,histsizeX+0.5,histsizeY+1,-0.5,histsizeY+0.5));
        //std::unique_ptr<TH3D> _h3(new TH3D("h3", "h3",histsizeX, 0, histsizeX, histsizeY,0, histsizeY, histsizeZ,0,histsizeZ));

        h1=std::move(_h1);
        h2=std::move(_h2);
        //  h3=std::move(_h3);

        auto now = std::chrono::system_clock::now();
        auto ms = std::chrono::time_point_cast<std::chrono::milliseconds>(now).time_since_epoch().count();
        Int_t seed = ms;
        std::unique_ptr<TRandom3> r(new TRandom3(seed));
        rndGen=std::move(r);

        std::cout<<"Start Analysis"<<std::endl;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------

TROOTAnalysis::~TROOTAnalysis(){
}


Int_t TROOTAnalysis::GetNofEntries(){
        return nofEntries;
}

void TROOTAnalysis::ApplyCut(Double_t c){

        for(Int_t i; i<Cevent->NHits(); i++) {

                if(Cevent->Hit(i)->EnergyDeposit() < c) {

                        Cevent->Hit(i)->SetEnergyDeposit(0.);


                }

        }

}


void TROOTAnalysis::PrintEdep(){

        for(Int_t i=0; i<Cevent->NHits(); i++) {
                std::cout<<Cevent->Hit(i)->EnergyDeposit()<<std::endl;
        }
}

// void TROOTAnalysis::PlotComparison(){
//
//   TFile * f1=new TFile("/home/iwsatlas1/emberger/Geant4/Data/ForwardCalor/CompletedAnalysis_BugFixed/ROOTFiles/_Pb1mm_absfirst_20x20mm_Res.root");
//
//   TCanvas * c1 =(TCanvas*)f1->Get("c1");
//   TGraph * g1 = new TGraph();
//   TList *l1=c1->GetListOfPrimitives();
//
//   g1=l1->Last();
//   TF1 *fc1=g1->GetFunction("rms");
//   Double_t standard_hundred=fc1->Eval(100);
//   Double_t standard_thousand= fc1->Eval(1000);
//
//
//
// }


void TROOTAnalysis::ResSummary(){
        //gStyle->SetOptTitle(kFALSE);
        //gStyle->SetPalette(kSolar);
//         TCanvas * rs1= new TCanvas("Angular Resolution Study", "Angular Resolution Study");
//         auto *mg= new TMultiGraph();
//         Double_t Tiledimension[5]={1,2,2.5,4,5};
//
//         // Double_t Pbabs_cut[4]={0.879, 0.909,1.104,1.507};
//         // Double_t pbgap_cut[4]={1.135, 1.201,1.382,1.687 };
//         // Double_t cuabs_cut[4]={0.463,0.352, 0.634,1.068 };
//         // Double_t cugap_cut[4]={0.745, 0.784, 0.969,1.303};
//         //
//         // Double_t noisePbabs_cut[4]={45.02,46.59,50.53,50.02};
//         // Double_t noisepbgap_cut[4]={46.08,47.18,50.42,49.83};
//         // Double_t noisecuabs_cut[4]={42.94,45.76,50.54,51.70};
//         // Double_t noisecugap_cut[4]={44.17,45.89,51.09,51.48};
//
//         Double_t Pbabs[5]={2.11, 2.17,2.22,2.375,2.464};
//         Double_t pbgap[5]={2.212, 2.277,2.318,2.456,2.526};
//         Double_t cuabs[5]={1.947,2.051,2.124 ,2.302,2.431 };
//         Double_t cugap[5]={2.05,2.121, 2.191,2.36,2.431 };
//
//         // Double_t noisePbabs[4]={45.24,46.3,50.13,49.25};
//         // Double_t noisepbgap[4]={46.1,47.15,50.37,49.22};
//         // Double_t noisecuabs[4]={42.93,45.8,51.19,52.04};
//         // Double_t noisecugap[4]={44.51,45.};
//
//
//
//
//         // for(Int_t i=0; i<4; i++) {
//         //         Pbabs[i]=Pbabs[i]-Pbabs_cut[i];
//         //         pbgap[i]=pbgap[i]-pbgap_cut[i];
//         //         cuabs[i]=cuabs[i]-cuabs_cut[i];
//         //         cugap[i]=cugap[i]-cugap_cut[i];
//         // }
//
//
//         TGraph * rsgPbAbs =  new TGraph(5, Tiledimension, Pbabs);
//         rsgPbAbs->SetTitle("Pb1mm_AbsorberFirst");
//         rsgPbAbs->GetXaxis()->SetTitle("dimension of scintillator tile [cm]");
//
//         rsgPbAbs->GetYaxis()->SetRangeUser(1,4);
//
//         rsgPbAbs->GetYaxis()->SetTitle("Angular Resolution Fit - stochastic term");
//
//         TGraph * rsgPbGap =  new TGraph(5,Tiledimension, pbgap);
//         rsgPbGap->SetTitle("Pb1mm_GapFirst");
//
//         TGraph * rsgCuAbs =  new TGraph(5, Tiledimension, cuabs);
//         rsgCuAbs->SetTitle("Cu2mm_AbsorberFirst");
//
//         TGraph * rsgCuGap =  new TGraph(5, Tiledimension, cugap);
//         rsgCuGap->SetTitle("Cu2mm_GapFirst");
//
//
//
//         rsgPbAbs->SetLineColor(kRed);
//         rsgPbAbs->SetLineStyle(1);
//         rsgPbAbs->SetMarkerStyle(23);
//         rsgPbAbs->SetLineWidth(2);
//         rsgPbAbs->SetMarkerSize(2);
//
//
//         rsgPbGap->SetLineColor(kRed);
//         rsgPbGap->SetLineStyle(10);
//         rsgPbGap->SetMarkerStyle(23);
//         rsgPbGap->SetMarkerSize(2);
//         rsgPbGap->SetLineWidth(3);
//
//         rsgCuAbs->SetLineColor(kBlue);
//         rsgCuAbs->SetLineStyle(1);
//         rsgCuAbs->SetMarkerStyle(23);
//         rsgCuAbs->SetLineWidth(2);
//         rsgCuAbs->SetMarkerSize(2);
//
//         rsgCuGap->SetLineColor(kBlue);
//         rsgCuGap->SetLineStyle(10);
//         rsgCuGap->SetMarkerStyle(23);
//         rsgCuGap->SetLineWidth(2);
//         rsgCuGap->SetMarkerSize(2);
//
//         //rs1->cd(0);
//         //
//         // mg->Add(rsgPbAbs, "");
//         // mg->Add(rsgPbGap, "");
//         // mg->Add(rsgCuAbs, "");
//         // mg->Add(rsgCuGap, "");
//
//         rsgPbAbs->Draw("APL");
//         rsgPbGap->Draw("PL");
//         rsgCuAbs->Draw("PL");
//         rsgCuGap->Draw("PL");
//
// // mg->Draw("");
//
// // gPad->BuildLegend();
//         auto * legend = new TLegend(0.1,0.7,0.48,0.9);
//         legend->AddEntry(rsgPbAbs, "1mmPb_AbsorberFirst", "l");
//         legend->AddEntry(rsgPbGap, "1mmPb_GapFirst", "l");
//         legend->AddEntry(rsgCuAbs, "2mmCu_AbsorberFirst", "l");
//         legend->AddEntry(rsgCuGap, "2mmCu_GapFirst", "l");
//
//         legend->Draw();
// TCanvas * rs2= new TCanvas("Angular_Resolution_Study", "Angular_Resolution_Study");
//
//         Double_t iX[2]={0.1, 1};
//
//
//         //Standard: 1mmLead_20x20mmTiles_AbsorberFirst
//         Double_t standard_hundred=0.2156;
//         Double_t standard_thousand=0.0636;
//
//         Double_t iYstandard[2]={1,1};
//
//
//         //CDR design is 1.75mmlead,AbsorberFirst
//         Double_t CDRPbabs_1cm[2]={ 0.2262/standard_hundred, 0.0688/standard_thousand};
//         Double_t CDRPbabs_2cm[2]={ 0.235/standard_hundred, 0.0719/standard_thousand};
//         Double_t CDRPbabs_25mm[2]={ 0.239/standard_hundred, 0.0737/standard_thousand};
//         Double_t CDRPbabs_4cm[2]={ 0.2565/standard_hundred, 0.0811/standard_thousand};
//         Double_t CDRPbabs_5cm[2]={ 0.2691/standard_hundred, 0.0853/standard_thousand};
//
//         Double_t hundredlayersPbabs2cm[2]={0.2173/standard_hundred,0.0632/standard_thousand};
//         Double_t hundredlayersCuabs2cm[2]={0.1993/standard_hundred,0.0573/standard_thousand};
//
//
//         Double_t iYPbabs_1cm[2]={ 0.2096/standard_hundred, 0.0618/standard_thousand};
//         Double_t iYPbgap_1cm[2]={ 0.22/standard_hundred, 0.0662/standard_thousand};
//         Double_t iYCuabs_1cm[2]={ 0.1928/standard_hundred, 0.0546/standard_thousand};
//         Double_t iYCugap_1cm[2]={ 0.2034/standard_hundred,  0.059/standard_thousand};
//
//
//         //Double_t iYPbabs_2cm[2]={ 0.2156/standard_hundred, 0.0636 /standard_thousand};
//         Double_t iYPbgap_2cm[2]={ 0.2266/standard_hundred, 0.0683/standard_thousand};
//         Double_t iYCuabs_2cm[2]={ 0.2032/standard_hundred, 0.0562/standard_thousand};
//         Double_t iYCugap_2cm[2]={ 0.2105/standard_hundred,  0.0612/standard_thousand};
//
//         Double_t iYPbabs_25mm[2]={ 0.2206/standard_hundred,  0.0651/standard_thousand};
//         Double_t iYPbgap_25mm[2]={ 0.231/standard_hundred, 0.0698/standard_thousand};
//         Double_t iYCuabs_25mm[2]={ 0.2104/standard_hundred, 0.0582/standard_thousand};
//         Double_t iYCugap_25mm[2]={ 0.2174/standard_hundred,  0.0628/standard_thousand};
//
//         Double_t iYPbabs_4cm[2]={ 0.2361/standard_hundred, 0.0701 /standard_thousand};
//         Double_t iYPbgap_4cm[2]={ 0.2445/standard_hundred, 0.0741/standard_thousand};
//         Double_t iYCuabs_4cm[2]={ 0.2281/standard_hundred, 0.0631/standard_thousand};
//         Double_t iYCugap_4cm[2]={ 0.234/standard_hundred,  0.0677/standard_thousand};
//
//         Double_t iYPbabs_5cm[2]={ 0.246/standard_hundred, 0.0768 /standard_thousand};
//         Double_t iYPbgap_5cm[2]={ 0.2525/standard_hundred, 0.0795/standard_thousand};
//         Double_t iYCuabs_5cm[2]={ 0.2414/standard_hundred, 0.0706/standard_thousand};
//         Double_t iYCugap_5cm[2]={ 0.246/standard_hundred, 0.0738/standard_thousand};
//
//
//
//
//         TGraph *gr1 =  new TGraph(2, iX, iYstandard);
//
//         gr1->GetYaxis()->SetRangeUser(0.8,1.1);
//         gr1->GetXaxis()->SetTitle("Energy [GeV]");
//         gr1->GetYaxis()->SetTitle("#frac{Resolution}{StandardResolution}");
//         gr1->GetYaxis()->SetTitleOffset(1.3);
//         gr1->SetTitle("Angular_Resolution_Study");
//
//         TGraph * gr3 =  new TGraph(2, iX, hundredlayersPbabs2cm);
//         //gr3->SetTitle("Pb1mm_AbsorberFirst");
//
//         TGraph * gr2 =  new TGraph(2,iX, iYCuabs_2cm);
//         //gr2->SetTitle("Pb1mm_GapFirst");
//
//         TGraph * gr5 =  new TGraph(2, iX, hundredlayersCuabs2cm);
//         //gr5->SetTitle("Cu2mm_AbsorberFirstiYCuabs_5cm
//
//         TGraph * gr4 =  new TGraph(2, iX, iYCugap_5cm);
//         //gr4->SetTitle("Cu2mm_GapFirst");
//
//
//
//         gr1->SetMarkerColor(kGray+2);
//         gr1->SetMarkerStyle(21);
//         gr1->SetMarkerSize(2);
//
//         gr2->SetMarkerColor(kOrange+7);
//         gr2->SetMarkerStyle(21);
//         gr2->SetMarkerSize(2);
//
//         gr3->SetMarkerColor(kGray+2);
//         gr3->SetMarkerStyle(20);
//         gr3->SetMarkerSize(2);
//
//         gr4->SetMarkerColor(kOrange+7);
//         gr4->SetMarkerStyle(20);
//         gr4->SetMarkerSize(2);
//
//         gr5->SetMarkerColor(kOrange+7);
//         gr5->SetMarkerStyle(20);
//         gr5->SetMarkerSize(2);
//
//         gr1->Draw("AP");
//         gr2->Draw("P");
//         gr3->Draw("P");
//         //gr4->Draw("P");
//         gr5->Draw("P");
//
//         auto * legend1 = new TLegend(0.1,0.7,0.48,0.9);
//         legend1->AddEntry(gr1, "1mmPb_AbsorberFirst_2cmTile_50layers", "p");
//
//         legend1->AddEntry(gr3, "1mmPb_AbsorberFirst_2cmTile_100layers", "p");
//         legend1->AddEntry(gr2, "2mmCu_AbsorberFirst_2cmTile_50layers", "p");
//         //legend1->AddEntry(gr4, "2mmCu_AbsorberFirst_2cmTile_100layers", "p");
//         legend1->AddEntry(gr5, "2mmCu_AbsorberFirst_2cmTile_100layers", "p");
//
//         legend1->Draw();

        // rs2->Print("/home/iwsatlas1/emberger/Geant4/Data/ForwardCalor/CompletedAnalysis_BugFixed/ComparisonPlots/100layerComparison.root");
        // rs2->Print("/home/iwsatlas1/emberger/Geant4/Data/ForwardCalor/CompletedAnalysis_BugFixed/ComparisonPlots/100layerComparison.C");
        // rs2->Print("/home/iwsatlas1/emberger/Geant4/Data/ForwardCalor/CompletedAnalysis_BugFixed/ComparisonPlots/100layerComparison.pdf");

        //
        // TCanvas * rs3 = new TCanvas();
        //
        // Double_t tilesize[5]={1,2,2.5,4,5};
        //
        // Double_t resolutionPbabs_hundred[5]={0.2096/standard_hundred,0.2156/standard_hundred,0.2206/standard_hundred,0.2361/standard_hundred,0.246/standard_hundred};
        // Double_t resolutionPbgap_hundred[5]={0.22/standard_hundred,0.2266/standard_hundred,0.231/standard_hundred,0.2445/standard_hundred,0.2525/standard_hundred};
        // Double_t resolutionCuabs_hundred[5]={0.1928/standard_hundred, 0.2032/standard_hundred,0.2104/standard_hundred,0.2281/standard_hundred,0.2414/standard_hundred};
        // Double_t resolutionCugap_hundred[5]={0.2034/standard_hundred,0.2105/standard_hundred,0.2174/standard_hundred,0.234/standard_hundred,0.246/standard_hundred};
        //   Double_t CDR_hundred[5]={ 0.2262/standard_hundred,0.235/standard_hundred,0.239/standard_hundred,0.2565/standard_hundred,0.2691/standard_hundred};
        //
        //
        // TGraph *gr6=new TGraph(5,tilesize,resolutionPbabs_hundred) ;
        // gr6->GetYaxis()->SetRangeUser(0.8,1.3);
        // gr6->GetXaxis()->SetTitle("Size of scintillator tile [cm]");
        // gr6->GetYaxis()->SetTitle("#frac{Resolution}{StandardResolution}");
        // gr6->GetYaxis()->SetTitleOffset(1.0);
        // gr6->SetTitle("Angular_Resolution_100MeV");
        //
        // TGraph * gr7=new TGraph(5,tilesize,resolutionPbgap_hundred) ;
        // TGraph * gr8=new TGraph(5,tilesize,resolutionCuabs_hundred) ;
        // TGraph * gr9=new TGraph(5,tilesize,resolutionCugap_hundred) ;
        // TGraph * gr10=new TGraph(5,tilesize, CDR_hundred);
        //
        // gr6->SetMarkerColor(kGray+2);
        // gr6->SetMarkerStyle(20);
        // gr6->SetMarkerSize(2);
        //
        // gr7->SetMarkerColor(kGray+2);
        // gr7->SetMarkerStyle(24);
        // gr7->SetMarkerSize(2);
        //
        // gr8->SetMarkerColor(kOrange+7);
        // gr8->SetMarkerStyle(20);
        // gr8->SetMarkerSize(2);
        //
        // gr9->SetMarkerColor(kOrange+7);
        // gr9->SetMarkerStyle(24);
        // gr9->SetMarkerSize(2);
        //
        // gr10->SetMarkerColor(kGray+2);
        // gr10->SetMarkerStyle(34);
        // gr10->SetMarkerSize(2);
        //
        // gr6->Draw("AP");
        // gr7->Draw("P");
        // gr8->Draw("P");
        // gr9->Draw("P");
        // gr10->Draw("P");
        //
        // auto * legend2 = new TLegend(0.1,0.7,0.48,0.9);
        // legend2->AddEntry(gr6, "1mmPb_AbsorberFirst", "p");
        // legend2->AddEntry(gr7, "1mmPb_ScintillatorFirst", "p");
        //         legend2->AddEntry(gr10, "1.75mmPb_AbsorberFirst_CDR", "p");
        // legend2->AddEntry(gr8, "2mmCu_AbsorberFirst", "p");
        // legend2->AddEntry(gr9, "2mmCu_ScintillatorFirst", "p");
        //
        //
        // legend2->Draw();
        //
        // rs3->Print("/home/iwsatlas1/emberger/Geant4/Data/ForwardCalor/CompletedAnalysis_BugFixed/ComparisonPlots/100MeVResolution_over_tilesize.root");
        // rs3->Print("/home/iwsatlas1/emberger/Geant4/Data/ForwardCalor/CompletedAnalysis_BugFixed/ComparisonPlots/100MeVResolution_over_tilesize.C");
        // rs3->Print("/home/iwsatlas1/emberger/Geant4/Data/ForwardCalor/CompletedAnalysis_BugFixed/ComparisonPlots/100MeVResolution_over_tilesize.pdf");


        // TCanvas * rs3 = new TCanvas();
        //
        // Double_t tilesize[5]={1,2,2.5,4,5};
        //
        // Double_t resolutionPbabs_thousand[5]={0.0618/standard_thousand,0.0636 /standard_thousand,0.0651/standard_thousand,0.0701 /standard_thousand,0.0768 /standard_thousand};
        // Double_t resolutionPbgap_thousand[5]={0.0662/standard_thousand,0.0683/standard_thousand,0.0698/standard_thousand,0.0741/standard_thousand,0.0795/standard_thousand};
        // Double_t resolutionCuabs_thousand[5]={0.0546/standard_thousand,0.0562/standard_thousand,0.0582/standard_thousand, 0.0631/standard_thousand,0.0706/standard_thousand};
        // Double_t resolutionCugap_thousand[5]={0.059/standard_thousand,0.0612/standard_thousand,0.0628/standard_thousand,0.0677/standard_thousand,0.0738/standard_thousand};
        // Double_t CDR_thousand[5]={0.0688/standard_thousand, 0.0719/standard_thousand, 0.0737/standard_thousand, 0.0811/standard_thousand, 0.0853/standard_thousand};
        //
        // TGraph *gr6=new TGraph(5,tilesize,resolutionPbabs_thousand) ;
        // gr6->GetYaxis()->SetRangeUser(0.8,1.4);
        // gr6->GetXaxis()->SetTitle("Size of scintillator tile [cm]");
        // gr6->GetYaxis()->SetTitle("#frac{Resolution}{StandardResolution}");
        // gr6->GetYaxis()->SetTitleOffset(1.3);
        // gr6->SetTitle("Angular_Resolution_1GeV");
        //
        // TGraph * gr7=new TGraph(5,tilesize,resolutionPbgap_thousand) ;
        // TGraph * gr8=new TGraph(5,tilesize,resolutionCuabs_thousand) ;
        // TGraph * gr9=new TGraph(5,tilesize,resolutionCugap_thousand) ;
        //  TGraph * gr10=new TGraph(5,tilesize, CDR_thousand);
        // gr6->SetMarkerColor(kGray+2);
        // gr6->SetMarkerStyle(20);
        // gr6->SetMarkerSize(2);
        //
        // gr7->SetMarkerColor(kGray+2);
        // gr7->SetMarkerStyle(24);
        // gr7->SetMarkerSize(2);
        //
        // gr8->SetMarkerColor(kOrange+7);
        // gr8->SetMarkerStyle(20);
        // gr8->SetMarkerSize(2);
        //
        // gr9->SetMarkerColor(kOrange+7);
        // gr9->SetMarkerStyle(24);
        // gr9->SetMarkerSize(2);
        //
        // gr10->SetMarkerColor(kGray+2);
        // gr10->SetMarkerStyle(34);
        // gr10->SetMarkerSize(2);
        //
        // gr6->Draw("AP");
        // gr7->Draw("P");
        // gr8->Draw("P");
        // gr9->Draw("P");
        // gr10->Draw("P");
        //
        // auto * legend2 = new TLegend(0.1,0.7,0.48,0.9);
        // legend2->AddEntry(gr6, "1mmPb_AbsorberFirst", "p");
        // legend2->AddEntry(gr7, "1mmPb_ScintillatorFirst", "p");
        //         legend2->AddEntry(gr10, "1.75mmPb_AbsorberFirst_CDR", "p");
        // legend2->AddEntry(gr8, "2mmCu_AbsorberFirst", "p");
        // legend2->AddEntry(gr9, "2mmCu_ScintillatorFirst", "p");
        //
        //
        // legend2->Draw();
        //
        // rs3->Print("/home/iwsatlas1/emberger/Geant4/Data/ForwardCalor/CompletedAnalysis_BugFixed/ComparisonPlots/1GeVResolution_over_tilesize.root");
        // rs3->Print("/home/iwsatlas1/emberger/Geant4/Data/ForwardCalor/CompletedAnalysis_BugFixed/ComparisonPlots/1GeVResolution_over_tilesize.C");
        // rs3->Print("/home/iwsatlas1/emberger/Geant4/Data/ForwardCalor/CompletedAnalysis_BugFixed/ComparisonPlots/1GeVResolution_over_tilesize.pdf");

}

void TROOTAnalysis::CalcAngularResolution(std::string pth){

        TCanvas * ca1=new TCanvas("c1","c1");

        std::string fileName;
        Int_t Estart=50;
        Double_t X[11];
        Double_t Y[11];

        for(Int_t i = 0; i<10; i++) {

                fileName=pth+"/gamma"+std::to_string(Estart)+"MeV_0.4:0.4_.root";
                std::cout<<fileName<<std::endl;
                TFile * f1=new TFile(fileName.c_str());

                TCanvas * c1 =(TCanvas*)f1->Get("ProjectionC");
                TH2D * hi1 = (TH2D*)c1->GetPrimitive("Prjection");

                TCanvas dummyC;

                Double_t rmsx=hi1->GetRMS(1);
                //Double_t rmsy=h1->GetRMS(2);
                X[i]=Estart;
                Y[i]=rmsx/1150;
                Estart+=100;

        }

        Estart=1500;
        fileName=pth+"/gamma"+std::to_string(Estart)+"MeV_0.4:0.4_.root";
        TFile * f1=new TFile(fileName.c_str());
        TCanvas * c1 =(TCanvas*)f1->Get("ProjectionC");
        TH2D * hi1 = (TH2D*)c1->GetPrimitive("Prjection");

        TCanvas dummyC;
        Double_t rmsx=hi1->GetRMS(1);
        //Double_t rmsy=h1->GetRMS(2);

        X[10]=Estart;
        Y[10]=rmsx/1150;

        TGraph * g1 = new TGraph(11, X, Y);
        std::string title=fileName.substr(61,23);
        g1->SetTitle(title.c_str());
        TF1 * rms_x = new TF1("rms", "sqrt(([0] * [0] / x) + ([2]*[2])  + ([1]*[1]/(x*x)))");
        //TF1 * rms_x = new TF1("rms", "[0]/(sqrt(x))");

        g1->Fit(rms_x);

        gStyle->SetOptFit(11111111);
        g1->GetXaxis()->SetTitle("Energy[MeV]");
        g1->GetYaxis()->SetTitle("#frac{RMS_x}{1000mm}");

        gStyle->SetOptFit();

        ca1->cd(0);
        g1->Draw("A*");

        std::string path=pth+"/"+title+"_Res.root";
        ca1->Print(path.c_str());

        path=pth+"/"+title+"_Res.C";
        ca1->Print(path.c_str());

        path=pth+"/"+title+"_Res.pdf";
        ca1->Print(path.c_str());


}

//void TROOTAnalysis::plotEvent(Int_t pev){     //plot 3DHisto of selected event
// TCanvas * plotcanvas1 = new TCanvas("eventplotter", "eventplotter");
//
// if(nofEntries<pev+1) {
//         std::cout<<"Tree has less than "<<pev+1<<" events, it has "<<nofEntries<<"."<<std::endl;
//         return;
// }
//
//
//
// //for(Int_t i=0;i<nent; i++){
//
// EcalTree->GetEntry(pev);
// Int_t pnh=Cevent->NHits();
//
// for(Int_t j=0; j<pnh; j++) {
//         h->Fill(Cevent->Hit(j)->X(), Cevent->Hit(j)->Y(),Cevent->Hit(j)->Z(), Cevent->Hit(j)->EnergyDeposit());
// }
//
// h->GetXaxis()->SetTitle("X");
// h->GetYaxis()->SetTitle("Y");
// h->GetZaxis()->SetTitle("Z");
//
// plotcanvas1->cd();
// h->Draw("BOX");
//
// std::string Path="/home/iwsatlas1/emberger/FuerFrank/PhotonEventDisplay/Event";
// std::string nr = std::to_string(pev);
//
// std::string extension = ".C";
// std::string PlotPath=Path+nr+extension;
// plotcanvas1->Print(PlotPath.c_str());
//
// extension = ".pdf";
// PlotPath=Path+nr+extension;
// plotcanvas1->Print(PlotPath.c_str());
// h->Reset();
//}

// void TROOTAnalysis::plotEventPion(Int_t pev){     //plot 3DHisto of selected event
//         TCanvas * plotcanvas1 = new TCanvas("eventplotter", "eventplotter");



//hA= new TH3D("ECalEvent1","ECalEvent1",histsizeX,0,histsizeX,histsizeY,0,histsizeY,histsizeZ,0,histsizeZ);
//hB= new TH3D("ECalEvent2","ECalEvent2",histsizeX,0,histsizeX,histsizeY,0,histsizeY,histsizeZ,0,histsizeZ);
//         if(nofEntries<pev+1) {
//                 std::cout<<"Tree has less than "<<pev+1<<" events, it has "<<nofEntries<<"."<<std::endl;
//                 return;
//         }
//
//
//
//         //for(Int_t i=0;i<nent; i++){
//
//         EcalTree->GetEntry(pev);
//         Int_t pnh=Cevent->NHits();
//
//         for(Int_t j=0; j<pnh; j++) {
//                 if(Cevent->Hit(j)->PhotNr()==1) {
//                         hA->Fill(Cevent->Hit(j)->X(), Cevent->Hit(j)->Y(),Cevent->Hit(j)->Z(), Cevent->Hit(j)->EnergyDeposit());
//                 }
//         }
//
//         //TColor * col1= new TColor(1,0,0,"red");
//
//         hA->GetXaxis()->SetTitle("X");
//         hA->GetYaxis()->SetTitle("Y");
//         hA->GetZaxis()->SetTitle("Z");
//         hA->SetMarkerColor(kRed);
//         hA->SetMarkerStyle(20);
//         hA->SetMarkerSize(0.3);
//         hA->Draw("");
//
//         for(Int_t j=0; j<pnh; j++) {
//                 if(Cevent->Hit(j)->PhotNr()==2) {
//                         hB->Fill(Cevent->Hit(j)->X(), Cevent->Hit(j)->Y(),Cevent->Hit(j)->Z(), Cevent->Hit(j)->EnergyDeposit());
//                 }
//         }
//
//         //TColor * col2 = new TColor(0,0,1,"blue");
//         hB->SetMarkerColor(kBlue);
//         hB->SetMarkerStyle(20);
//         hB->SetMarkerSize(0.3);
//         hB->GetXaxis()->SetTitle("X");
//         hB->GetYaxis()->SetTitle("Y");
//         hB->GetZaxis()->SetTitle("Z");
//         hB->Draw("SAME");
//
//         std::string Path="/home/iwsatlas1/emberger/PionPlots/Event";
//         std::string nr = std::to_string(pev);
//
//         std::string extension = ".C";
//         std::string PlotPath=Path+nr+extension;
//         plotcanvas1->Print(PlotPath.c_str());
//
//         // extension = ".png";
//         // PlotPath=Path+nr+extension;
//         // plotcanvas1->Print(PlotPath.c_str());
//
//         extension = ".pdf";
//         PlotPath=Path+nr+extension;
//         plotcanvas1->Print(PlotPath.c_str());
//
//         hB->Reset();
//         hA->Reset();
// }


//-----------------------------------------------------------------------------------------------------------------------------------------------------


void TROOTAnalysis::PlotRMSx(){
        TCanvas * rms=new TCanvas("rmsx_over_E");


        // Double_t y[10]={0.5298,0.3507,0.3624,0.2524,0.222,0.1861,0.228,0.1751,0.1425,0.1481}; // 0.0 slope
        // Double_t x[10]={0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1};
        //
        Double_t y[10]={0.5398,0.3428,0.3087,0.2319,0.2117,0.1933,0.2243,0.1716,0.1411,0.1588};                                                                      //0.4slope
        Double_t x[10]={0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1};
        TF1 * rms_x = new TF1("rms", "sqrt(([0] * [0] / x) + ([2]*[2])  + ([1]*[1]/(x*x)))");

        TGraph *gRMS1=new TGraph(10, x, y );
        gRMS1->Fit(rms_x);
        gStyle->SetOptFit(11111111);
        gRMS1->GetXaxis()->SetTitle("Energy[GeV]");
        gRMS1->GetYaxis()->SetTitle("#frac{RMS_x}{100mm}");
        gRMS1->GetYaxis()->SetTitleOffset(1.3);
        rms->cd(0);
        gRMS1->Draw("A*");
}


//---------------------------------------------------------------------------------------------------------------------------------------------------------

// void TROOTAnalysis::SampleFromHIst(){
//
//         TCanvas * dummy = new TCanvas();
//         dummy->Divide(2,2,0.01,0.01);
//         TFile *f = new TFile("/home/iwsatlas1/emberger/Geant4/Current/4piSensitiveDetector_mod/B4-build/pi0_kinematics.root");
//         //f.ls();
//         Double_t energy;
//         Double_t cos_theta;
//
//         TH2D * hist1 = (TH2D*)f->Get("pi0_energyAngle");
//
//         TH1D * sample1= new TH1D("sample1", "sample1", 100, 0, 1);
//         TH1D * sample2= new TH1D("sample2", "sample2", 200, -1, 1);
//
//         TRandom3 * rndGen = new TRandom3();
//
//         TVector3 BeamDirection(0.,0.,1.);
//         BeamDirection=BeamDirection.Unit();
//
//         TH1D * check = new TH1D("Check", "check", 200, -1, 1);
//         TH2D * fd = new TH2D("fd", "fd", 200,-2,2,200,-2,2);
//         for(Int_t i=0; i<1000000; i++) {
//                 hist1->GetRandom2(energy, cos_theta);
//                 sample1->Fill(energy);
//                 sample2->Fill(cos_theta);
//
//
//                 if(cos_theta>0) {
//                         Double_t theta = TMath::ACos(cos_theta);
//
//                         Double_t radius = TMath::Tan(theta);
//
//                         Double_t a=rndGen->Uniform(0,TMath::TwoPi());
//
//                         Double_t x=TMath::Cos(a)*TMath::Abs(radius);
//                         Double_t y=TMath::Sin(a)*TMath::Abs(radius);
//                         fd->Fill(x,y);
//                         TVector3 MomentumDirection(x+BeamDirection.x(), y+BeamDirection.y(), BeamDirection.z());
//                         MomentumDirection=MomentumDirection.Unit();
//                         check->Fill(TMath::Cos(MomentumDirection.Angle(BeamDirection)));
//                 }
//                 else if(cos_theta==0) {
//                         //Double_t theta = TMath::ACos(cos_theta);
//
//                         Double_t a=rndGen->Uniform(0,TMath::TwoPi());
//
//                         Double_t x=TMath::Cos(a)*1;
//                         Double_t y=TMath::Sin(a)*1;
//                         fd->Fill(x,y);
//
//                         TVector3 MomentumDirection(x+BeamDirection.x(), y+BeamDirection.y(), 0);
//                         MomentumDirection=MomentumDirection.Unit();
//                         check->Fill(TMath::Cos(MomentumDirection.Angle(BeamDirection)));
//                 }
//                 else if(cos_theta<0) {
//                         Double_t theta = TMath::ACos(cos_theta);
//
//                         Double_t radius = TMath::Tan(theta);
//
//                         Double_t a=rndGen->Uniform(0,TMath::TwoPi());
//
//                         Double_t x=TMath::Cos(a)*TMath::Abs(radius);
//                         Double_t y=TMath::Sin(a)*TMath::Abs(radius);
//                         fd->Fill(x,y);
//
//                         TVector3 MomentumDirection(x+BeamDirection.x(), y+BeamDirection.y(), -BeamDirection.z());
//                         MomentumDirection=MomentumDirection.Unit();
//                         check->Fill(TMath::Cos(MomentumDirection.Angle(BeamDirection)));
//                 }
//
//         }
//
//         dummy->cd(1);
//         sample1->Draw();
//         dummy->cd(2);
//         sample2->Draw();
//         dummy->cd(3);
//         check->Draw();
//         dummy->cd(4);
//         fd->Draw();
//
//
// }



void TROOTAnalysis::PrintERes(){

        TCanvas * res = new TCanvas("EnergyResolution");
        //TCanvas * gap = new TCanvas("GapEnergy", "GapEnergy");

        //Double_t y[12]={3.945/21.74,6.73/64.4,8.89/107.2,10.56/150,12.18/193,13.38/235.8,14.5/278.4,15.96/321.4, 16.57/364.1,17.32/407.2,22.07/642.2 };
        Double_t x[2]={0.1,1.0};

        Double_t y_pbabs[2]={0.108, 0.03};
        Double_t y_pbgap[2]={0.131, 0.033};
        Double_t y_cuabs[2]={0.115,0.035};
        Double_t y_cugap[2]={0.148,0.041 };
        Double_t y_cdr[2]={0.16, 0.0516};

        TGraph * re1 = new TGraph(2, x, y_pbabs);
        re1->GetYaxis()->SetRangeUser(0.02,0.18);
        re1->GetXaxis()->SetTitle("Energy [GeV]");
        re1->GetYaxis()->SetTitle("Energy Resolution");
        re1->GetYaxis()->SetTitleOffset(1.3);
        re1->SetTitle("EnergyResolution");

        TGraph * re2= new TGraph(2,x,y_pbgap);
        TGraph * re3= new TGraph(2,x,y_cuabs);
        TGraph * re4= new TGraph(2,x,y_cugap);
        TGraph * re5= new TGraph(2,x,y_cdr);

        re1->SetMarkerColor(kGray+2);
        re1->SetMarkerStyle(20);
        re1->SetMarkerSize(2);

        re2->SetMarkerColor(kGray+2);
        re2->SetMarkerStyle(24);
        re2->SetMarkerSize(2);

        re3->SetMarkerColor(kOrange+7);
        re3->SetMarkerStyle(20);
        re3->SetMarkerSize(2);

        re4->SetMarkerColor(kOrange+7);
        re4->SetMarkerStyle(24);
        re4->SetMarkerSize(2);

        re5->SetMarkerColor(kGray+2);
        re5->SetMarkerStyle(34);
        re5->SetMarkerSize(2);

        re1->Draw("AP");
        re2->Draw("P");
        re3->Draw("P");
        re4->Draw("P");
        re5->Draw("P");

        auto * legend2 = new TLegend(0.1,0.7,0.48,0.9);
        legend2->AddEntry(re1, "1mmPb_AbsorberFirst", "p");
        legend2->AddEntry(re2, "1mmPb_ScintillatorFirst", "p");
        legend2->AddEntry(re5, "1.75mmPb_AbsorberFirst_CDR", "p");
        legend2->AddEntry(re3, "2mmCu_AbsorberFirst", "p");
        legend2->AddEntry(re4, "2mmCu_ScintillatorFirst", "p");

        legend2->Draw();

        res->Print("/home/iwsatlas1/emberger/Geant4/Data/ForwardCalor/CompletedAnalysis_BugFixed/Plots/ERes_100layers/EnergyResolutionComparison.root");
        res->Print("/home/iwsatlas1/emberger/Geant4/Data/ForwardCalor/CompletedAnalysis_BugFixed/Plots/ERes_100layers/EnergyResolutionComparison.C");
        res->Print("/home/iwsatlas1/emberger/Geant4/Data/ForwardCalor/CompletedAnalysis_BugFixed/Plots/ERes_100layers/EnergyResolutionComparison.pdf");


        // TImage * img1 = TImage::Create();
        // img1->FromPad(res);
        // //img1->Scale(1000, 1000);
        // img1->WriteImage("/home/iwsatlas1/emberger/Geant4/Current/SensitiveDetector/B4-build/B4c/GammaEnergyandSlopeScan_Analysis/ERes.png");
        // delete img1;
        //
        // Double_t Egun[12]={50,100,200,300,400,500,700,1000,2500,5000,7500,10000};
        // Double_t Egap[12]={19.25,38.47,77.51,115.9,154.2,192.9,268.8,384,956.2,1905,2836,3779};
        //
        // Double_t errx[12]={0,0,0,0,0,0,0,0,0,0,0,0};
        // Double_t erry[12]={4.55, 5.942, 9.056, 12.06, 11.59, 13.71, 17.98, 20.25, 53.25, 69.23, 94.6, 146.7};
        //
        // TF1 * Efficiency = new TF1("Efficiency", "[0]*x+[1]");
        //
        // TGraphErrors * gap1 = new TGraphErrors(12, Egun, Egap, errx, erry );
        // gap1->SetTitle("Efficiency");
        // gap1->Fit(Efficiency);
        // gap1->GetXaxis()->SetTitle("Gun Energy[MeV]");
        // gap1->GetYaxis()->SetTitle("Gap Energy[MeV]");
        // gap1->GetYaxis()->SetTitleOffset(1.3);
        // gStyle->SetOptFit();
        //
        // gap->cd(0);
        // gap1->Draw("A*");
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------


void TROOTAnalysis::PlotProjection(Double_t distance, Int_t event){
        if(!flg) {

                Double_t k=(-distance-EstimatePhoton1[event].first.Z())/EstimatePhoton1[event].second.Z();

                TVector3 Intersect(EstimatePhoton1[event].first.X()+k*EstimatePhoton1[event].second.X(),
                                   EstimatePhoton1[event].first.Y()+k*EstimatePhoton1[event].second.Y(),
                                   EstimatePhoton1[event].first.Z()+k*EstimatePhoton1[event].second.Z());


                projection_pca->Fill(Intersect.X(), Intersect.Y());
                
        }

        else {
                Double_t k=(-distance-DirectionPhoton1[event].first.Z())/DirectionPhoton1[event].second.Z();

                TVector3 Intersect(DirectionPhoton1[event].first.X()+k*DirectionPhoton1[event].second.X(),
                                   DirectionPhoton1[event].first.Y()+k*DirectionPhoton1[event].second.Y(),
                                   DirectionPhoton1[event].first.Z()+k*DirectionPhoton1[event].second.Z());

                projection_minimization->Fill(Intersect.X(), Intersect.Y());

        }

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------
Bool_t TROOTAnalysis::PCAEvent(Int_t event){

        Int_t nHits = Cevent->NHits();

        std::vector<TVector3> eigenVecs;

        std::vector<Double_t> showerCOG(3, 0.);

        //calculate the average shower position

        // TVector3 dummyvec(0.,0.,0.);
        // std::vector<TVector3> sample;
        // nHits=1000;
        // for(Int_t i=0;i<nHits; i++) {
        //         dummyvec.SetX(dummyvec.X()+0.1);
        //         dummyvec.SetY(dummyvec.Y()+rndGen->Gaus(0.0, 0.2));
        //         dummyvec.SetZ(dummyvec.Z()+rndGen->Gaus(0.0, 0.2));
        //
        //         sample.push_back(dummyvec);
        // }

        Double_t eges=0.;
        for(Int_t i=0; i<nHits; i++) {

                auto hitnormal=TransformCoordinates(Cevent->Hit(i)->X(), Cevent->Hit(i)->Y(), Cevent->Hit(i)->Z());
                Double_t edep=Cevent->Hit(i)->EnergyDeposit();

                std::tuple<Double_t, Double_t, Double_t> hit=std::make_tuple(std::get<2>(hitnormal), std::get<1>(hitnormal), std::get<0>(hitnormal));


                showerCOG.at(0)+=(std::get<0>(hit)*edep);
                showerCOG.at(1)+=(std::get<1>(hit)*edep);
                showerCOG.at(2)+=(std::get<2>(hit)*edep);
                eges+=edep;

                // showerCOG.at(0)+=sample[i].X();
                // showerCOG.at(1)+=sample[i].Y();
                // showerCOG.at(2)+=sample[i].Z();

        }

        showerCOG.at(0)/=eges;
        showerCOG.at(1)/=eges;
        showerCOG.at(2)/=eges;

        // showerCOG.at(0)/=nHits;
        // showerCOG.at(1)/=nHits;
        // showerCOG.at(2)/=nHits;
        //
        // showerCOG.at(0)/=1000;
        // showerCOG.at(1)/=1000;
        // showerCOG.at(2)/=1000;


        Double_t xi2 = 0.;
        Double_t xiyi = 0.;
        Double_t xizi = 0.;
        Double_t yi2 = 0.;
        Double_t yizi = 0.;
        Double_t zi2 = 0.;
        Double_t weightSum = 0.;

        //Get symmetric matrix


        for(Int_t i=0; i<nHits; i++) {

                auto hitnormal=TransformCoordinates(Cevent->Hit(i)->X(), Cevent->Hit(i)->Y(), Cevent->Hit(i)->Z());

                std::tuple<Double_t, Double_t, Double_t> hit=std::make_tuple(std::get<2>(hitnormal), std::get<1>(hitnormal), std::get<0>(hitnormal));

                Double_t weight=Cevent->Hit(i)->EnergyDeposit();
                //Double_t weight= 1.;

                Double_t x=(std::get<0>(hit)-showerCOG.at(0))*weight;
                Double_t y=(std::get<1>(hit)-showerCOG.at(1))*weight;
                Double_t z=(std::get<2>(hit)-showerCOG.at(2))*weight;

                // Double_t x=(sample[i].X()-showerCOG.at(0))*weight;
                // Double_t y=(sample[i].Y()-showerCOG.at(1))*weight;
                // Double_t z=(sample[i].Z()-showerCOG.at(2))*weight;

                weightSum+=weight * weight;
                xi2 += x * x;
                xiyi += x * y;
                xizi += x * z;
                yi2 += y * y;
                yizi += y * z;
                zi2 += z * z;
        }

        Eigen::Matrix3f mat;

        mat <<  xi2,    xiyi,   xizi,
                xiyi,   yi2,    yizi,
                xizi,   yizi,   zi2;

        mat /= weightSum;


        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> decomp(mat);

        if (decomp.info() == Eigen::ComputationInfo::Success) {
                using eigenValColPair = std::pair<float, size_t>;
                std::vector<eigenValColPair> eigenValColVec = {
                        eigenValColPair(decomp.eigenvalues()(0), 0),
                        eigenValColPair(decomp.eigenvalues()(1), 1),
                        eigenValColPair(decomp.eigenvalues()(2), 2)
                };

                std::sort(eigenValColVec.begin(), eigenValColVec.end(),[](const eigenValColPair& left, const eigenValColPair& right){
                        return left.first > right.first;
                });

                Eigen::Matrix3f eigenVectors(decomp.eigenvectors());
                for (const auto &pair : eigenValColVec) {
                        TVector3 temp(
                                eigenVectors(0, pair.second),
                                eigenVectors(1, pair.second),
                                eigenVectors(2, pair.second));

                        eigenVecs.push_back(temp);
                }
                EstimatePhoton1.push_back(std::make_pair(TVector3(showerCOG.at(0), showerCOG.at(1), showerCOG.at(2)), eigenVecs[0]));
                return true;
        }
        else {
                std::cerr << "WARNING: PCA decompose failure in event "<<event <<"."<<"Number of hits: "<<nHits<<std::endl;
                return false;
        }



}
//--------------------------------------------------------------------------------------------------------------

void TROOTAnalysis::CalcCOGPion(Int_t event){            //calculate vector of (X,Y) tuples containing layerwise center of gravity



        // Variables for fit
        Double_t xerr=0;
        Double_t yerr=0;
        Double_t cgx=0;
        Double_t cgy=0;
        Double_t cgz=0;
        Double_t Eweight=0;

        // Int_t segment;
        //
        Double_t EGesPH1=0;
        // Double_t EGesPH2=0;
        EcalTree->GetEntry(event);         //grab event from tree
        Eges = Cevent->GapEnergy();
        Int_t cnh = Cevent->NHits();

        Double_t integral;
        Int_t nofLayerswithEdep=0;

        for(Int_t i=0; i<nofLayers; i++) {         //loop over all layers in event

                for(Int_t j =0; j<cnh; j++) {         //loop over all hits in laver i

                        //std::cout<<Cevent->Hit(j)->PhotNr()<<std::endl;
                        Double_t edep= Cevent->Hit(j)->EnergyDeposit();
                        if(Cevent->Hit(j)->Z()==i &&  edep > 0.0) {

                                //std::cout<<"fill: "<<Cevent->Hit(j)->PhotNr()<<"from photon mode:"<<phnr<<std::endl;
                                h1->Fill(Cevent->Hit(j)->X(), Cevent->Hit(j)->Y(), edep);

                        }         //fill histogram with hits of layer i
                }
                //std::cout<<photonNR<<" : "<<std::endl;
                integral=h1->Integral();

                if(integral!=0) {         //check if layer containes energy
                        nofLayerswithEdep++;

                        Eweight=integral/Cevent->GapEnergy();  //calculate weight

                        EGesPH1+=integral;

                        cgx=h1->GetMean(1);
                        xerr=h1->GetMeanError(1);
                        //extract center of gravity and error
                        cgy=h1->GetMean(2);
                        yerr=h1->GetMeanError(2);

                        if(yerr<0.00001) {yerr=1/TMath::Sqrt(12); }
                        if(xerr<0.00001) {xerr=1/TMath::Sqrt(12); }

                        cgz=i;

                        //std::cout<<"segment: "<<std::endl;

                        std::tuple<Double_t, Double_t, Double_t> COGtransformed;

                        //h3->Fill(cgx,cgy,cgz);

                        COGtransformed=TransformCoordinates(cgx, cgy, cgz /*, segment*/);

                        auto cg=std::make_tuple(std::get<0>(COGtransformed), std::get<1>(COGtransformed), std::get<2>(COGtransformed), xerr, yerr, Eweight);

                        //  std::cout<<"------------------------"<<std::endl;

                        coglist.push_back(cg);


                }         // end of clustering


                //ClusteredHits.clear();

                h1->Reset();

                h2->Reset();

        }                 //end of event



        //showerCOGPhoton1.push_back(TransformCoordinates(h3->GetMean(1), h3->GetMean(2), h3->GetMean(3) /*, segment*/));
        COGCollectionPH1.push_back(coglist);
        EnergyPhoton1.push_back(std::make_pair(EGesPH1/0.52, 1 /*segment*/));
        coglist.clear();

        //h3->Reset();
        // showerdepthH->Fill(nofLayerswithEdep);


// totalEnergyH->Fill((EGesPH1+EGesPH2)/1000);

//check->cd(0);
//h3->GetYaxis()->SetRange(40,60);
//h3->Draw("*");
}

//-------------------------------------------------------------------------------------------------------------------------------


std::tuple<Double_t, Double_t, Double_t > TROOTAnalysis::TransformCoordinates(Double_t x, Double_t y, Double_t z /*, Double_t seg*/){

        std::tuple<Double_t,Double_t, Double_t> transformedCOGs;


        // if(seg==1) {
        Double_t offsetX=-1*calsizeXY/2+tiledimX/2;
        Double_t offsetY=-1*calsizeXY/2+tiledimY/2;
        Double_t offsetZ;
        if(absofirst) {
                offsetZ=AbsoThickness+ GapThickness/2;
                //std::cout<<"true"<<std::endl;
        }
        else{
                offsetZ=GapThickness/2;
                //  std::cout<<"false"<<std::endl;

        }
        //  std::cout<<"X: "<<offsetX<<"Y:"<<offsetY<<"Z: "<<offsetZ<<std::endl;

        TVector3 coordinates((offsetX + x * tiledimX),                      // X
                             (offsetY + y * tiledimY),                      // Y
                             (offsetZ + z * (AbsoThickness+GapThickness))); // Z

        //coordinates.Print();

        //rotate into correct segment

        transformedCOGs=std::make_tuple(coordinates.X(), coordinates.Y(), coordinates.Z());


        return transformedCOGs;
}


//-------------------------------------------------------------------------------------------------------------------------------------


void TROOTAnalysis::FitCOGsPion(Int_t event){
        gErrorIgnoreLevel = 1001;
        EcalTree->GetEntry(event);

        TVector3 ph1_orig=Cevent->MomentumPh1().Unit();

        Fcn myfcn;
        myfcn.SetMode("photon3D");

        MnUserParameters upar;
        double error_minimizer_parameters = 1e-4;


        myfcn.SetCOGs(COGCollectionPH1[event]);

        upar.Add("a_1", /*EstimatePhoton1[event].first.X()*/ 0., error_minimizer_parameters);
        upar.Add("a_2", /*EstimatePhoton1[event].first.Y()*/ 0., error_minimizer_parameters);
        upar.Add("a_3", /*EstimatePhoton1[event].first.Z()*/ 0., error_minimizer_parameters);

        // upar.Fix("a_1");
        // upar.Fix("a_2");
        // upar.Fix("a_3");

        //std::cout<<"set COGPH1"<<std::endl;
        upar.Add("v_1", /*EstimatePhoton1[event].second.X()*/ 0., error_minimizer_parameters);
        upar.Add("v_2", /*EstimatePhoton1[event].second.Y()*/ 0., error_minimizer_parameters);
        upar.Add("v_3", /*EstimatePhoton1[event].second.Z()*/ 1., error_minimizer_parameters);

        //  upar.SetLimits("v_3", 0.6, 1.0);
        //cout << upar << endl;

        MnMigrad migrad(myfcn, upar, 2);

        myfcn.SetCurrentEvent(event);

        //std::cout<<"event: "<<events<<std::endl;
        FunctionMinimum min = migrad();

        //std::cout<<min<<std::endl;

        MnUserParameterState uParState = min.UserState();

        //std::cout<<"attempt adding params for photon 1"<<std::endl;

        TVector3 A_1(uParState.Value("a_1"), uParState.Value("a_2"), uParState.Value("a_3"));
        TVector3 V_1(uParState.Value("v_1"), uParState.Value("v_2"), uParState.Value("v_3"));
        V_1=V_1.Unit();

        auto tp1 = std::make_pair(A_1, V_1);
        //std::cout<<"errga"<<std::endl;
        DirectionPhoton1.push_back(tp1);
        //std::cout<<"adding params for photon 1"<<std::endl;

        TVector3 true_direction=TVector3(0.4,0.4,1.0).Unit();

        dx1->Fill(DirectionPhoton1[event].second.X());
        dy1->Fill(DirectionPhoton1[event].second.Y());
        dz1->Fill(TMath::Sqrt(DirectionPhoton1[event].second.Z()*DirectionPhoton1[event].second.Z()));

}

void TROOTAnalysis::DrawHists(){
//

        dx1->GetXaxis()->SetTitle("X component of Direction Vector");
        dx1->GetYaxis()->SetTitle("Entries");

        dx2->GetXaxis()->SetTitle("X component of Direction Vector");
        dx2->GetYaxis()->SetTitle("Entries");

        dy1->GetXaxis()->SetTitle("Y component of Direction Vector");
        dy1->GetYaxis()->SetTitle("Entries");

        dy2->GetXaxis()->SetTitle("Y component of Direction Vector");
        dy2->GetYaxis()->SetTitle("Entries");

        dz1->GetXaxis()->SetTitle("Z component of Direction Vector");
        dz1->GetYaxis()->SetTitle("Entries");

        dz2->GetXaxis()->SetTitle("Z component of Direction Vector");
        dz2->GetYaxis()->SetTitle("Entries");

        D->Divide(3,2,0.01,0.01);
        gStyle->SetOptStat(2222);
        D->cd(1);
        dx1->Draw();
        D->cd(2);
        dy1->Draw();
        D->cd(3);
        dz1->Draw();
        D->cd(4);
        dx2->Draw();
        D->cd(5);
        dy2->Draw();
        D->cd(6);
        dz2->Draw();

        projectionC->Divide(3,1,0.01,0.01);

        projectionC->cd(1);
        projection_pca->Draw("colz");

        projectionC->cd(2);
        projection_minimization->Draw("colz");

        projectionC->cd(3);
        projection_correlation->Draw("colz");

        std::string path=savepath+"/"+filename+"_.root";
        projectionC->Print(path.c_str());

        // path=savepath+"/"+filename+".pdf";
        // projectionC->Print(path.c_str());

        // showerdepth->cd(0);
        // showerdepthH->Draw();
//
//         Delta1->Divide(2,2,0.01,0.01);
//
//         gStyle->SetOptStat(2222);
//
//         delx->GetXaxis()->SetTitle("DeltaX[mm]");
//         dely->GetXaxis()->SetTitle("DeltaY[mm]");
//         delz->GetXaxis()->SetTitle("DeltaZ[mm]");
//
//         delx->GetYaxis()->SetTitle("Entries");
//         dely->GetYaxis()->SetTitle("Entries");
//         delz->GetYaxis()->SetTitle("Entries");
//         appdist1->GetXaxis()->SetTitle("Distance of closest Approach [mm]");
//         appdist1->GetYaxis()->SetTitle("Entries");
//
//         Delta1->cd(1);
//         delx->Draw();
//         Delta1->cd(2);
//         dely->Draw();
//         Delta1->cd(3);
//         delz->Draw();
//         Delta1->cd(4);
//         appdist1->Draw();
//
//         // Delta2->cd(0);
//         // delx->Draw();
//         //
//         // Delta3->cd(0);
//         // appdist1->Draw();
//         //
//         // Delta4->cd(0);
//         // delz->Draw();
//
//
//         // histograms for different invariaont mass calculations
//
//         InvMassReco->GetXaxis()->SetTitle("InvariantMass[MeV]");
//         InvMassReco->GetYaxis()->SetTitle("Entries");
//
//         InvMassSim->GetXaxis()->SetTitle("InvariantMass[MeV]");
//         InvMassSim->GetYaxis()->SetTitle("Entries");
//
//         InvMassAllreco->GetXaxis()->SetTitle("InvariantMass[MeV]");
//         InvMassAllreco->GetYaxis()->SetTitle("Entries");
//
//         //histograms for Angle_vs_InvariantMass plots
//
//         Angle_vs_massReco->GetXaxis()->SetTitle("Invariant Mass [MeV]");
//         Angle_vs_massReco->GetYaxis()->SetTitle("Angle[deg]");
//         Angle_vs_massReco->GetYaxis()->SetTitleOffset(1.4);
//
//         Angle_vs_massSim->GetXaxis()->SetTitle("Invariant Mass [MeV]");
//         Angle_vs_massSim->GetYaxis()->SetTitle("Angle[deg]");
//         Angle_vs_massSim->GetYaxis()->SetTitleOffset(1.4);
//
//
//         Angle_vs_massAllreco->GetXaxis()->SetTitle("Invariant Mass [MeV]");
//         Angle_vs_massAllreco->GetYaxis()->SetTitle("Angle[deg]");
//         Angle_vs_massAllreco->GetYaxis()->SetTitleOffset(1.4);
//
//         IM1->Divide(2,2,0.01,0.01);
//
//         gStyle->SetOptStat(2222);
//         // IM1->cd(1);
//         // InvMassReco->Draw();
//
//         IM1->cd(1);
//         InvMassSim->Draw();
//
//         IM1->cd(2);
//         InvMassAllreco->Draw();
//
//
//
//         Int_t palette[5];
//         palette[0] = 3;
//         palette[1] = 4;
//         palette[2] = 5;
//         palette[3] = 6;
//         palette[4] = 2;
//         gStyle->SetPalette(5, palette);
//         gStyle->SetOptStat(0);
//
//         // IM1->cd(4);
//         // Angle_vs_massReco->GetZaxis()->SetNdivisions(5);
//         // Angle_vs_massReco->Draw("colz");
//         // Angle_vs_massReco->SetMaximum(5);
//         // Angle_vs_massReco->SetMinimum(0);
//
//         IM1->cd(3);
//         Angle_vs_massSim->GetZaxis()->SetNdivisions(5);
//         Angle_vs_massSim->Draw("colz");
//         Angle_vs_massSim->SetMaximum(5);
//         Angle_vs_massSim->SetMinimum(0);
//
//         IM1->cd(4);
//         Angle_vs_massAllreco->GetZaxis()->SetNdivisions(5);
//         Angle_vs_massAllreco->Draw("colz");
//         Angle_vs_massAllreco->SetMaximum(5);
//         Angle_vs_massAllreco->SetMinimum(0);
//         gStyle->SetOptStat(2222);
//         PiLoc->Divide(3,3,0.01,0.01);
//
//         pilocX->GetXaxis()->SetTitle("Closest Approach Position X [mm]");
//         pilocX->GetYaxis()->SetTitle("Entries");
//         PiLoc->cd(4);
//         pilocX->Draw();
//
//         pilocY->GetXaxis()->SetTitle("Closest Approach Position Y [mm]");
//         pilocY->GetYaxis()->SetTitle("Entries");
//         PiLoc->cd(5);
//         pilocY->Draw();
//
//         pilocZ->GetXaxis()->SetTitle("Closest Approach Position Z [mm]");
//         pilocZ->GetYaxis()->SetTitle("Entries");
//         PiLoc->cd(6);
//         pilocZ->Draw();
//
//         inputX->GetXaxis()->SetTitle("Closest Approach Position X [mm]");
//         inputX->GetYaxis()->SetTitle("Entries");
//         PiLoc->cd(1);
//         inputX->Draw();
//
//         inputY->GetXaxis()->SetTitle("Closest Approach Position Y [mm]");
//         inputY->GetYaxis()->SetTitle("Entries");
//         PiLoc->cd(2);
//         inputY->Draw();
//
//         inputZ->GetXaxis()->SetTitle("Closest Approach Position Z [mm]");
//         inputZ->GetYaxis()->SetTitle("Entries");
//
//         PiLoc->cd(3);
//         inputZ->Draw();
//
//         // PiLoc1->cd(0);
//         // inputZ->Draw();
//         //
//         // PiLoc2->cd(0);
//         // pilocZ->Draw();
//
//         PiLoc->cd(7);
//         invmasFitX->Draw();
//
//         PiLoc->cd(8);
//         invmasFitY->Draw();
//
//         PiLoc->cd(9);
//         invmasFitZ->Draw();
//
//         delR->GetXaxis()->SetTitle("Distance of true to reconstructed vertex [mm]");
//         delR->GetYaxis()->SetTitle("Entries");
//
//         DeltaR->cd(0);
//         delR->Draw();
//         Double_t x, q;
//         q = 0.5; // 0.5 for "median"
//         delR->ComputeIntegral(); // just a precaution
//         delR->GetQuantiles(1, &x, &q);
//         std::cout << "median = " << x << std::endl;
//         DeltaR->Update();
//
//         TLine *line = new TLine(x,0,x,16);
//         line->SetLineColor(kRed);
//         line->SetLineWidth(3);
//         line->Draw();
//
//         DeltaR_fit->cd(0);
//         delR_fit->Draw();
//         delR_fit->ComputeIntegral(); // just a precaution
//         delR_fit->GetQuantiles(1, &x, &q);
//         std::cout << "median = " << x << std::endl;
//         DeltaR_fit->Update();
//
//         TLine *line_fit = new TLine(x,0,x,16);
//         line_fit->SetLineColor(kRed);
//         line_fit->SetLineWidth(3);
//         line_fit->Draw();
//
//         std::string path="/home/iwsatlas1/emberger/Presentations/CALICE_japan/4Pi/X_700mm_400MeV/";
//
//         std::string path1=path+"DeltasDirection.pdf";
//         D->Print(path1.c_str());
//         path1=path+"DeltasDirection.C";
//         D->Print(path1.c_str());
//
//         std::string path2=path+"Invariant_mass.pdf";
//         IM1->Print(path2.c_str());
//         path2=path+"Invariant_mass.C";
//         IM1->Print(path2.c_str());
//
//         std::string path3=path+"Pi0Vertex.pdf";
//         Delta1->Print(path3.c_str());
//         path3=path+"Pi0Vertex.C";
//         Delta1->Print(path3.c_str());
//
//         std::string path4=path+"PionLocator.pdf";
//         PiLoc->Print(path4.c_str());
//         path4=path+"PionLocator.C";
//         PiLoc->Print(path4.c_str());
//
//
//         // path="/home/iwsatlas1/emberger/Presentations/CALICE_japan/Finals/";
//         //
//         // std::string path5=path+"X_VertexREco.pdf";
//         // Delta2->Print(path5.c_str());
//         // std::string path6=path+"Dist_VertexREco.pdf";
//         // Delta3->Print(path6.c_str());
//         // std::string path7=path+"Z_VertexREco.pdf";
//         // Delta4->Print(path7.c_str());
//         //
//         // std::string path8=path+"PionLocator_inputZ.pdf";
//         // PiLoc1->Print(path8.c_str());
//         // std::string path9=path+"PionLocatorZ.pdf";
//         // PiLoc2->Print(path9.c_str());
//         //
//         // std::string path10=path+"3DDist.pdf";
//         // DeltaR->Print(path10.c_str());
//         //
//         // std::string path11=path+"3DDist_fit.pdf";
//         // DeltaR_fit->Print(path11.c_str());
//
//
//         ang1->Divide(1,2,0.01,0.01);
//         ang1->cd(1);
//         angh_true->Draw();
//
//         ang1->cd(2);
//         angh_rec->Draw();
//
//         totalEnergyC->cd(0);
//         totalEnergyH->GetXaxis()->SetTitle("GeV");
//         totalEnergyH->GetYaxis()->SetTitle("Entries");
//
//         totalEnergyH->Draw();
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------

//void TROOTAnalysis::PlotChimap(Int_t event){

// TCanvas * ChiMap=new TCanvas("ChiMapC","Chimap");
// ChiMap->Divide(2,1,0.01,0.01);
// TH3D * chimap = new TH3D("chimap", "chimap",1000,-500,500,1000,-500,500,1000,-1600,400);
//
//
//
// for(Int_t i=0; i<nofEntries; i++) {
//
//         EcalTree->GetEntry(i);
//
//         // Double_t e1=Cevent->EnergyPhoton1();
//         // Double_t e2=Cevent->EnergyPhoton2();
//         Double_t e1=EnergyPhoton1[i];
//         Double_t e2=EnergyPhoton2[i];
//         TVector3 startvalues1(std::get<0>(ClosestApproach[i]) - std::get<0>(showerCOGPhoton1[i]),
//                               std::get<1>(ClosestApproach[i]) - std::get<1>(showerCOGPhoton1[i]),
//                               std::get<2>(ClosestApproach[i]) - std::get<2>(showerCOGPhoton1[i]));
//         Double_t deltaT1=startvalues1.Mag();
//
//         TVector3 startvalues2(std::get<0>(ClosestApproach[i]) - std::get<0>(showerCOGPhoton2[i]),
//                               std::get<1>(ClosestApproach[i]) - std::get<1>(showerCOGPhoton2[i]),
//                               std::get<2>(ClosestApproach[i]) - std::get<2>(showerCOGPhoton2[i]));
//         Double_t deltaT2=startvalues2.Mag();
//
//         TVector3 PointPh1(std::get<0>(FitParamsPions[i][0]) * std::get<2>(ClosestApproach[i]) + std::get<1>(FitParamsPions[i][0]),
//                           std::get<2>(FitParamsPions[i][0]) * std::get<2>(ClosestApproach[i]) + std::get<3>(FitParamsPions[i][0]),
//                           std::get<2>(ClosestApproach[i]));
//
//         TVector3 PointPh2(std::get<0>(FitParamsPions[i][1]) * std::get<2>(ClosestApproach[i]) + std::get<1>(FitParamsPions[i][1]),
//                           std::get<2>(FitParamsPions[i][1]) * std::get<2>(ClosestApproach[i]) + std::get<3>(FitParamsPions[i][1]),
//                           std::get<2>(ClosestApproach[i]));
//
//         PointPh1.Print();
//         PointPh2.Print();
//         std::cout<<std::get<0>(ClosestApproach[i])<<":"<<std::get<1>(ClosestApproach[i])<<":"<<std::get<2>(ClosestApproach[i])<<std::endl;
//
//         Double_t chisq;
//
//         for(Int_t x = -500; x<501; x+=2) {
//                 for(Int_t y = -500; y<501; y+=2) {
//                         for(Int_t z = -1600; z< 401; z+=2) {
//
//                                 Double_t sigPh1=(1/TMath::Sqrt(e1))*deltaT1;
//                                 Double_t sigPh2=(1/TMath::Sqrt(e2))*deltaT2;
//
//
//                                 Double_t delR1=TMath::Sqrt(((PointPh1.X()-x)*(PointPh1.X()-x)) + ((PointPh1.Y()-y)*(PointPh1.Y()-y)) + ((PointPh1.Z()-z)*(PointPh1.Z()-z)));
//                                 Double_t delR2=TMath::Sqrt(((PointPh2.X()-x)*(PointPh2.X()-x)) + ((PointPh2.Y()-y)*(PointPh2.Y()-y)) + ((PointPh2.Z()-z)*(PointPh2.Z()-z)));
//
//                                 TVector3 Dir_ph1(x-std::get<0>(showerCOGPhoton1[i]),y-std::get<1>(showerCOGPhoton1[i]),z-std::get<2>(showerCOGPhoton1[i]));
//
//                                 Dir_ph1=Dir_ph1.Unit();
//
//                                 TLorentzVector Lv_ph1(Dir_ph1.X()*e1,Dir_ph1.Y()*e1,Dir_ph1.Z()*e1,e1);
//
//
//                                 TVector3 Dir_ph2(x-std::get<0>(showerCOGPhoton2[i]),y-std::get<1>(showerCOGPhoton2[i]),z-std::get<2>(showerCOGPhoton2[i]));
//
//                                 Dir_ph2=Dir_ph2.Unit();
//
//                                 TLorentzVector Lv_ph2(Dir_ph2.X()*e2,Dir_ph2.Y()*e2,Dir_ph2.Z()*e2,e2);
//
//                                 TLorentzVector Lv_pi=Lv_ph1+Lv_ph2;
//
//                                 Double_t Curr_invmass=Lv_pi.M();
//
//
//
//
//                                 chisq=((delR1/sigPh1)*(delR1/sigPh1))+((delR2/sigPh2)*(delR2/sigPh2))+(((Curr_invmass-134.9766)/5)*((Curr_invmass-134.9766)/5));
//
//                                 chimap->Fill(x,y,z,chisq);
//
//
//                         }
//                         //std::cout<<y<<std::endl;
//                 }
//
//         }
//         std::cout<<"projection"<<std::endl;
//         ChiMap->cd(1);
//         chimap->Project3D("xz")->Draw("colz");
//         std::cout<<"projection"<<std::endl;
//         ChiMap->cd(2);
//         chimap->Project3D("yz")->Draw("colz");
//
//         std::string ChiMapPath="Chi2Mapz/ChiMap" +std::to_string(Cevent->GapEnergy())+"_" +std::to_string(Cevent->GunPos().Z())+"_"+ std::to_string(i) + ".pdf";
//
//         ChiMap->Print(ChiMapPath.c_str());
//
//         ChiMapPath="Chi2Mapz/ChiMap" +std::to_string(Cevent->GapEnergy()) +"_" +std::to_string(Cevent->GunPos().Z())+"_"+ std::to_string(i) + ".C";
//
//         ChiMap->Print(ChiMapPath.c_str());
//
//         chimap->Reset();
//
//
//
//
// }

//}
