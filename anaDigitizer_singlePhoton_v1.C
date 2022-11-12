#include "anaDigitizer_update.h"
// Program to calculate amplitude and charge from the
// TODO

void anaDigitizer_singlePhoton_v1(UInt_t prcEvents = 1000){

  const UInt_t eventSz=1024; // Size of each event
  //****************************************************************
  // Signal to analize inicialization
  //****************************************************************
  float  buffer[eventSz]; //Buffer to store the event
  FILE *ptr; //Pointer to the file to process
  const Char_t * file = "/Users/solangelmac/cernbox/Students/Janek/2300VSmallPMTatt5/wave_0.dat"; // File to analize
  // Parameters to analize the signal
  Double_t Wstart =95; // Start of the time expected of the signal. Begining of the time window to analize the signal.
  Double_t Wend = 130; // End of the time expected of the signal. End of the time window to analize the signal.
  Double_t thrs = 10; // Leading edge threshold
  Int_t sign = -1 ;// Sign fo the signal to analize

  // Number of entries for the analized file
  Int_t dataPointSz = sizeof buffer[0];
  cout <<  "Size of data point: " << sizeof buffer[0] <<" Bytes" << endl;
  Int_t fSize = filesize(file);
  UInt_t nEvents = fSize/(eventSz*dataPointSz);
  //nEvents = 50;
  printf("File size %i Bytes\nEvent size (%i)\nTotal Events(%i): %s\n\n", fSize, eventSz*dataPointSz, nEvents, file);
  ptr = fopen(file,"rb");  // r for read, b for binary


  //****************************************************************
  // TRIGGER 0 signal inicialization
  //****************************************************************
  float  bufferTrg0[eventSz]; //Buffer to store the event
  FILE *ptrTrg0; //Pointer to the file to process
  const Char_t * fileTrg0 = "//Users/solangelmac/cernbox/Students/Janek/2300VSmallPMTatt5/TR_0_0.dat"; // Trigger: time reference signal
  // Parameters to analize the signal
  Double_t WstartTrg0 =10; // Start of the time expected of the signal. Begining of the time window to analize the signal.
  Double_t WendTrg0 = 400; // End of the time expected of the signal. End of the time window to analize the signal.
  Double_t thrsTrg0 = 6; // Leading edge threshold
  Int_t signTrg0 = -1 ;// Sign fo the signal to analize

  // Number of entries for the time reference file
  dataPointSz = sizeof bufferTrg0[0];
  cout <<  "Size of data point (Trg0): " << sizeof bufferTrg0[0] <<" Bytes" << endl;
  Int_t fSizeTrg0 = filesize(fileTrg0);
  UInt_t nEventsTrg0 = fSizeTrg0/(eventSz*dataPointSz);
  nEvents = nEventsTrg0;
  if(prcEvents>0)  nEvents = prcEvents;

  printf("File size (Trg0) %i Bytes\nEvent size (%i)\nTotal Events(%i): %s\n\n", fSizeTrg0, eventSz*dataPointSz, nEventsTrg0, fileTrg0);
  ptrTrg0 = fopen(fileTrg0,"rb");  // r for read, b for binary


  //****************************************************************
  //****************************************************************

  //Fill in nseconds the time buffer for all the signals
  Double_t dt[eventSz];
  for (UInt_t i = 0; i < eventSz; i++) {dt[i]=i*timeRes;}

  Double_t amplitude;
  float RawAmp;

  // Plots for the signals
  TGraph* grAv  = nullptr;// For the analized signal
  TGraph* grRaw = new TGraph(eventSz);
  TGraph* grAvTrg0 = nullptr;// For the trigger 0 signal
  TGraph* grRawTrg0 = new TGraph(eventSz);

  // Filling histograms and plots
  for(UInt_t event=0; event < nEvents; event++){
    // Plots for the signals
    grAv = new TGraph(eventSz);
    grAv -> GetYaxis() -> SetRange(-1000,1000);
    grAvTrg0 = new TGraph(eventSz);
    grAvTrg0 -> GetYaxis() -> SetRange(-1000,1000);

    // Put one signal (1024*UInt32) in the buffer
    fread(buffer,sizeof(buffer),1,ptr);
    fread(bufferTrg0,sizeof(bufferTrg0),1,ptrTrg0);

    // Producing the graphs
    // The TRIGGER 0 first to calculate the time shift from the other signal
    h1Type h1ParsTrg0 = GetBaseLine(bufferTrg0,50);
    for (UInt_t i = 0; i < eventSz; i++){
      grAvTrg0 ->SetPoint(i,dt[i],bufferTrg0[i]*ampRes-h1ParsTrg0.mean);
    }
    grType sgPropTrg0 = GetGrProp(grAvTrg0, thrsTrg0, signTrg0, WstartTrg0, WendTrg0);


    // Graph from the signal input
    h1Type h1Pars = GetBaseLine(buffer,50);
    for (UInt_t i = 0; i < eventSz; i++){
      amplitude = buffer[i]*ampRes - h1Pars.mean;
      // Fill histo and plot with time shifted
      grAv ->SetPoint(i,dt[i]-sgPropTrg0.zcross,amplitude);
      h2SignalShift -> Fill(dt[i]-sgPropTrg0.zcross,amplitude);
      // Not shifted time histogram
      h2Signal -> Fill(dt[i],amplitude);
    }//

    // Get signal properties
    grType sgProp = GetGrProp(grAv, thrs, sign, Wstart, Wend);
    h1ampSig->Fill(sign*sgProp.vmin);
    h1charge->Fill(GetCharge(grAv,105,30));
    // Quality control histograms
    h1BaseLineAll->Fill(h1Pars.mean);
    h1RMSAll->Fill(h1Pars.rms);

    //keep the last graph for debuging
    if(event < nEvents-1){
      delete grAv;
      delete grAvTrg0;
      }
    else{
      //Filling the amplitude raw values to debug and quality control
      for (UInt_t i = 0; i < eventSz; i++){
        RawAmp = buffer[i];//&0x0000fff;
        //printf("point %d, buffer %d\n", i,buffer[i]);
        grRaw ->SetPoint(i,dt[i],RawAmp);
        grAv ->SetPoint(i,dt[i],RawAmp*ampRes- h1Pars.mean);
        grRawTrg0 ->SetPoint(i,dt[i],bufferTrg0[i]);
        grAvTrg0 ->SetPoint(i,dt[i],bufferTrg0[i]*ampRes-h1ParsTrg0.mean);
      }
    }// Keep last hsto for debuging END

  } // Fill histograms and plots

  TCanvas * c1 = new TCanvas("c1","c1 title",1600,800);
  c1->Divide(3,1);
  c1->cd(1);
  grAv ->SetLineColor(kBlack);
  grAv ->SetLineWidth(1);
  grAv ->Draw("ALP*");

  //***************************************************************************************************
  //***************************************************************************************************
  // CODE FOR DEBUGGING PRUPOSSES
  grType sgPropTrg0 = GetGrProp(grAvTrg0, thrsTrg0, signTrg0, WstartTrg0, WendTrg0);
  //grType sgPropTrg1 = GetGrProp(grAvTrg1, thrsTrg1, signTrg1, WstartTrg1, WendTrg1);
    h1Type h1Pars = GetBaseLine(buffer,50);
    for (UInt_t i = 0; i < eventSz; i++){
      amplitude = buffer[i]*ampRes - h1Pars.mean;
      // Fill histo and plot with time shifted
      grAv ->SetPoint(i,dt[i]-sgPropTrg0.zcross,amplitude);
      h2SignalShift -> Fill(dt[i]-sgPropTrg0.zcross,amplitude);
      // Not shifted time histogram
      h2Signal -> Fill(dt[i],amplitude);
    }//
  grType sgProp = GetGrProp(grAv, thrs, sign, Wstart, Wend);  // Graph from the signal input
  //***************************************************************************************************
  //***************************************************************************************************

  c1->cd(2);
  grAvTrg0 ->SetLineColor(kBlue);
  grAvTrg0 ->SetLineWidth(1);
  grAvTrg0 -> SetMarkerColor(kBlue);
  grAvTrg0 ->Draw("ALP*");
  c1->cd(3);

  TCanvas * cRaw = new TCanvas("cRaw","cRaw title",1400,800);
  grRaw ->SetLineColor(kBlack);
  grRaw ->SetLineWidth(1);
  grRaw ->Draw("ALP*");

  TCanvas * c2 = new TCanvas("c2","c2 title",1400,800);
  c2 -> Divide(2,1);
  c2 -> cd(1) -> SetLogz();
  h2Signal -> Draw("COLZ");
  c2 -> cd(2) ->SetLogz();
  h2SignalShift -> Draw("COLZ");

  TCanvas * c3 = new TCanvas("c3","Base line for all events",1400,800);
  c3 -> Divide (2,1);
  c3 -> cd(1) -> SetLogy();
  h1BaseLineAll -> Draw();
  c3 -> cd(2) -> SetLogy();
  h1RMSAll -> Draw();

  TCanvas * c4 = new TCanvas("c4","charge",800,800);
  c4-> SetLogy();
  h1charge->Draw();

  TCanvas * c5 = new TCanvas("c5","Amplitude",800,800);
  c5 -> Divide (2,1);
  c5 -> cd(1) -> SetLogy();
  c5 -> cd(2) -> SetLogy();
  h1ampSig->Draw();
}
