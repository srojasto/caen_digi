#include "TF1.h"
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include <fstream>
#include "TCanvas.h"

#include <iostream>
#include <cstdio>

// histograms

TH2 * h2SignalShift = new TH2D("h2SignalShift", "Signal Average shifted; time (ns);Amplitude (mV);",1024*3,-0.2*1024,0.2*1023*2 ,4096*0.24414,-4095*0.24414,4096*0.24414);
TH2 * h2Signal = new TH2D("h2Signal", "Signal Average; time (ns);Amplitude (mV);", 1024,0,0.2*1023 ,4096*2,-4096*0.24414,4095*0.24414);

TH1 * h1BaseLineAll = new TH1D("h1BaseLineAll", "Base Line from all events;Base line level (mV);Entries",1024,0,1024);
TH1 * h1RMSAll = new TH1D("h1RMSAll", "RMS from all events;RMS (mV);Entries",1000,0,100);
TH1 * h1charge = new TH1D("h1charge", ";Charge (fC);Entries",11000,-100,1000);
TH1 * h1ampSig = new TH1D("h1ampSig", "Signal;Amplitude (mV);Entries",2000,-1000,1000);
TH1 * h1ampTr1 = new TH1D("h1ampTr1", "Trigger;Amplitude (mV);Entries",2000,-1000,1000);

// *************************************************
//        Variables and structure declarations
// *************************************************
Double_t ampRes = 0.244140625;// 1/4096;// amplitude resolution = 0.2 mV per bit (LSB) ;
Double_t timeRes = 0.2;// 1/5e9 samples/second -> 0.2 nsecons per bit ;

struct h1Type {
    Double_t mean;
    Double_t rms;
};

// Structure to extract the information of each event
struct grType {
  UInt_t n;
  Double_t * y;// = gr->GetY();
  Double_t * x;// = gr->GetY();
  UInt_t locmax;// = TMath::LocMax(n,y);
  UInt_t locmin;// = TMath::LocMax(n,y);
  Double_t vmax;// = y[locmax];
  Double_t vmin;// = y[locmax];
  Double_t tmin;// = x[locmax];
  Double_t tmax;// = x[locmax];
  Double_t tshift;// = y[locmax];
  UInt_t loctshift;
  Bool_t UnderThr = kFALSE;
  Double_t zcross; // Zero crossing
};


// *************************************************
//        Function prototype
// *************************************************
std::ifstream::pos_type filesize(const char* filename);
h1Type GetBaseLine(float * tempBuffer);
grType GetGrProp(TGraph * gr, Double_t threshold, Int_t PulseSign , Double_t WStart, Double_t WEnd );


// *************************************************
//        Functions
// *************************************************
std::ifstream::pos_type filesize(const char* filename)
{
    std::ifstream in(filename, std::ifstream::ate | std::ifstream::binary);
    return in.tellg();
}


h1Type GetBaseLine(float * tempBuffer, const UInt_t baseLineSample=100){
  TH1D * h1BaseLine = new TH1D("h1BaseLine", "Base Line;Base line level (mV);Entries",1024,0,1024);
  for (UInt_t i = 0; i < baseLineSample; i++){ //Calculate base line
      h1BaseLine -> Fill(tempBuffer[i]*ampRes);
  }

  h1Type h1Result;
  h1Result.mean = h1BaseLine -> GetMean();
  h1Result.rms = h1BaseLine -> GetRMS();
  delete h1BaseLine;
  return h1Result;
}

grType GetGrProp(TGraph * gr, Double_t threshold = 0, Int_t PulseSign = -1, Double_t WStart =0, Double_t WEnd = 204.8){ // Graph properties
  grType graph;

  // Assign vertical (x) and horizontal (y) of the graph to pointers X and Y
  graph.n = gr ->GetN();
  graph.y = gr ->GetY();
  graph.x = gr ->GetX();

  // Locate array points for the maximum and minimum values
  graph.locmax = TMath::LocMax(graph.n,graph.y);
  graph.locmin = TMath::LocMin(graph.n,graph.y);

  // Get value of the minimum and maximum voltages and time
  graph.vmax = graph.y[graph.locmax];
  graph.vmin = graph.y[graph.locmin];
  graph.tmin = graph.x[graph.locmin];
  graph.tmax = graph.x[graph.locmax];

  // Finding the location of a point over specific threshold
  for(UInt_t i = 0; i < graph.n; i++){
    if( (graph.x[i] > WStart) && (graph.x[i] < WEnd) ){
      // The fit is automatized based on a treshold and a range to analize the signal
      // Fit a pol1 funtion to get the zero crossing of a positive signals
      if(graph.y[i] >= threshold && PulseSign == 1){
        // Fit a pol1 function to the signals
        TF1 *f = new TF1("f", "[0]+[1]*x", graph.x[i], graph.x[i]+3);
        gr->Fit("f","R");

        graph.zcross = -f->GetParameter(0)/f->GetParameter(1);
        //cout<< "Trigger   zero crossing time = "<< graph.zcross << endl;
        break;
      }
      //if(graph.y[i] <= threshold && PulseSign == -1){
      if(PulseSign == -1){
        // Fit a pol1 function to the signals
        // TF1 *f = new TF1("f", "[0]+[1]*x", graph.x[i], graph.x[i]+3);
        // gr->Fit("f","R");
        //cout<< "Tmin "<< graph.tmin<< endl;
        TF1 *f = new TF1("f", "[0]+[1]*x", graph.tmin-8,graph.tmin);
        gr->Fit("f","RQ");
        graph.zcross = -f->GetParameter(0)/f->GetParameter(1);
        //cout<< "Trigger   zero crossing time = "<< graph.zcross << endl;
        break;
      }

      if(graph.y[i] < threshold && PulseSign == -1){
        graph.tshift = graph.x[i];
        graph.loctshift = i;
        //printf("time %f, %f, %f\r", graph.x[i], graph.y[i], PulseSign * threshold);
        break;
      }

      if(graph.y[i] > threshold && PulseSign == 1){
        graph.tshift = graph.x[i];
        graph.loctshift = i;
        //printf("time %f, %f, %f \r", graph.x[i], graph.y[i], PulseSign * threshold);
        break;
      }
    }
    if(graph.x[i] > WStart && graph.y[i] < (PulseSign * threshold)) graph.UnderThr=kTRUE;

    // TODO handle signals below the threshold
  }
  return graph;
}


double GetCharge(TGraph * gr, double  WStart=0 ,   double  WEnd = 204.8)
{
  // Integrate in the interval
  Double_t charge=0;
  Double_t amplitude=0;
  for(int i = 0; i < gr ->GetN(); i++){
    amplitude = -gr->GetPointY(i);
    if( (gr->GetPointX(i)>WStart) && (gr->GetPointX(i)<WEnd) ){
      charge += amplitude*(0.2/50);
    }
  }
  //cout<<"Charge "<< charge<< endl;
  return charge;
}
