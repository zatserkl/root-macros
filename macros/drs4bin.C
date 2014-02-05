/*
g++ -Wall drs4bin.cpp
./a.out drs4bin.cpp
*/

#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>
#include <TCanvas.h>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cstdarg>
#include <cstdlib>
#include <cstdio>

using std::cout;      using std::endl;

/*
TODO: 
1) exect EHDR and C001,C002, etc.
2) read complete 4-bytes word and union with two 2-bytes words
3) exceptions for readout error
4) easy way to set brach addresses to read the tree (MakeClass may work)
   -- how to do with different number of channels?
*/

class DRS4bin {
public:
   bool status;
   bool operator !() {return !status;}
   static Float_t t[1024];
   static Float_t v1[1024];
   static Float_t v2[1024];
   static Float_t v3[1024];
   static Float_t v4[1024];
   static Bool_t usedchan[4];
   static Int_t event;
   static Int_t year;
   static Int_t month;
   static Int_t day;
   static Int_t hour;
   static Int_t minute;
   static Int_t second;
   static Int_t millisecond;
   std::string ifname;
   std::ifstream ifile;
   Long64_t ifsize;
   TFile* ofile;
   TTree* tree;

   static TGraph* plot(Int_t ch)
   {
      if (!usedchan[ch]) return 0;
      const Int_t color[] = {2, 4, 6, 8};

      Float_t* v = 0;
      if (ch == 0) v = v1;
      if (ch == 1) v = v2;
      if (ch == 2) v = v3;
      if (ch == 3) v = v4;
      TGraph* gr = new TGraph(1024,t,v);
      gr->SetNameTitle(Form("gr_evt_%d_chan_%d",event,channel), Form("gr_evt_%d_chan_%d",event,channel));
      gr->SetMarkerStyle(6);
      gr->SetMarkerColor(color[ch]);
      gr->Draw("ap");
      return gr;
   }

   DRS4bin(std::string ifname_0): status(false)
                                  , ifname(ifname_0)
                                  , ofile(0)
                                  , tree(0)
   {
      ifile.open(ifname.c_str(), std::ios::binary);
      if (!ifile) {
         cout<< "File not found: " << ifname <<endl;
         return;
      }
      cout<< "processing file " << ifname <<endl;

      //-- file length
      // TODO 
      ifile.seekg(0, std::ios::end);
      ifsize = ifile.tellg();
      ifile.seekg(0);
      if (ifsize < 0) {
         cout<< "input file error: ifsize = " << ifsize <<endl;
         ifile.close();
         return;
      }
      // else cout<< "ifsize = " << ifsize <<endl;

      status = true;
      usedchan[0] = usedchan[1] = usedchan[2] = usedchan[3] = kFALSE;
   }

   //-- Convert
   void Convert(Int_t entry1=0, Int_t entry2=0)
   {
      // set exception mask for ifile
      ifile.exceptions(std::ifstream::failbit | std::ifstream::badbit);

      // read the first event to figure out used channels
      try {
         readDRS4event();
      }
      catch (std::ifstream::failure e) {
         //cout<< "read error: " << e.what() <<endl;
         cout<< "DRS4bin::Convert: exception caught: read error at the first event" <<endl;
         return;
      }

      cout<< "DRS4 channels: ";
      for (int ich=0; ich<4; ++ich) cout<< ich+1 << (usedchan[ich]? ": yes ": ": no ");
      cout<<endl;

      ofile = TFile::Open(Form("%s.root",ifname.c_str()),"recreate");
      tree = new TTree("drs", "produced from DRS4 binary file");
      tree->SetMarkerStyle(6);
      tree->SetMarkerColor(46);
      tree->SetLineColor(46);

      tree->Branch("t", &t, "t[1024]/F");
      // book used channels
      if (usedchan[0]) tree->Branch("v1", &v1, "v1[1024]/F");
      if (usedchan[1]) tree->Branch("v2", &v2, "v2[1024]/F");
      if (usedchan[2]) tree->Branch("v3", &v3, "v3[1024]/F");
      if (usedchan[3]) tree->Branch("v4", &v4, "v4[1024]/F");
      // rest of the event
      tree->Branch("usedchan", &usedchan, "usedchan[4]/b");
      tree->Branch("event", &event, "event/I");
      tree->Branch("year", &year, "year/I");
      tree->Branch("month", &month, "month/I");
      tree->Branch("day", &day, "day/I");
      tree->Branch("hour", &hour, "hour/I");
      tree->Branch("minute", &minute, "minute/I");
      tree->Branch("second", &second, "second/I");
      tree->Branch("millisecond", &millisecond, "millisecond/I");

      if (entry1 == 0) tree->Fill();

      // TODO calculate event size and seekg to proper position

      // process the rest of the data
      try {
         Int_t ientry = 0;
         bool res = true;
         while (res) {
            res = readDRS4event();
            if (res) {
               ++ientry;
               if (ientry >= entry1) tree->Fill();
               if (entry2 > 0 && ientry >= entry2) break;
            }
         }
      }
      catch (std::ifstream::failure e) {
         //cout<< "read error: " << e.what() <<endl;
         cout<< "DRS4bin::Convert: exception caught: read error at event " << tree->GetEntries() <<endl;
      }

      cout<< "write " << tree->GetEntries() << " entries into file " << ofile->GetName() <<endl;
      ofile->Write();

      new TCanvas;
      tree->Draw("v1:t","");
   }

   //-- readDRS4event()
   bool readDRS4event()
   {
      Long64_t pos = ifile.tellg();
      if (pos < 0) return false;

      if (pos == ifsize) return false;

      union {
         Char_t buffer[24];
         struct {
            char title[4];
            UInt_t event;        // serial number starting from 1
            UShort_t year;
            UShort_t month;
            UShort_t day;
            UShort_t hour;
            UShort_t minute;
            UShort_t second;
            UShort_t millisecond;
            UShort_t reserved;
         };
      } header;

      ifile.read(header.buffer, sizeof(header.buffer));

      // sanity check
      if (false
         or (header.title[0] != 'E')
         or (header.title[1] != 'H')
         or (header.title[2] != 'D')
         or (header.title[3] != 'R')
         )
      {
         // this is not an Event Header
         ifile.seekg(pos);
         // cout<< "this is not an Event Header: ";
         // cout.write(header.title,4);
         // cout<<endl;
         return false;
      }

      event = header.event;
      year = header.year;
      month = header.month;
      day = header.day;
      hour = header.hour;
      minute = header.minute;
      second = header.second;
      millisecond = header.millisecond;

      //cout<< "year " << year << " month " << month << " day " << day << " hour " << hour << " minute " << minute << " second " << second << " millisecond " << millisecond <<endl;

      // read time bins
      char* t_buffer = (char*) t;
      ifile.read(t_buffer, 1024*4);

      // read channels
      bool res = true;
      while (res) {
         res = readDRS4channel();
      }

      return true;
   }

   //-- readDRS4channel
   bool readDRS4channel()
   {
      Long64_t pos = ifile.tellg();
      if (pos < 0) return false;

      if (pos == ifsize) {return false;}

      char header[4];

      ifile.read(header, sizeof(header));

      Float_t* v = 0;
      int channel;
      switch (header[3]) {
         case '1': channel = 1; v = v1; usedchan[0] = kTRUE; break;
         case '2': channel = 2; v = v2; usedchan[1] = kTRUE; break;
         case '3': channel = 3; v = v3; usedchan[2] = kTRUE; break;
         case '4': channel = 4; v = v4; usedchan[3] = kTRUE; break;
      }

      if (v == 0 or (header[0] != 'C') or (header[1] != '0') or (header[2] != '0'))
      {
         // this is not a channel header
         ifile.seekg(pos);    // return to previous position
         // cout<< "this is not a Channel Header: ";
         // cout.write(header,4);
         // cout<<endl;
         return false;
      }

      union {
         Char_t byte[2048];
         UShort_t word2[1024];
      } voltage;

      ifile.read(voltage.byte, sizeof(voltage.byte));
      for (int ipoint=0; ipoint<1024; ++ipoint) {
         v[ipoint] = (Float_t(voltage.word2[ipoint]) - 32767.5) / 65536.;
      }

      return true;
   }
};

// define static members
Float_t DRS4bin::t[1024];
Float_t DRS4bin::v1[1024];
Float_t DRS4bin::v2[1024];
Float_t DRS4bin::v3[1024];
Float_t DRS4bin::v4[1024];
Bool_t DRS4bin::usedchan[4];
Int_t DRS4bin::event;
Int_t DRS4bin::year;
Int_t DRS4bin::month;
Int_t DRS4bin::day;
Int_t DRS4bin::hour;
Int_t DRS4bin::minute;
Int_t DRS4bin::second;
Int_t DRS4bin::millisecond;

//void drs4bin(const char* ifname="data_meander.bin", Int_t entry1=0, Int_t entry2=0)
void drs4bin(const char* ifname, Int_t entry1=0, Int_t entry2=0)
{
   DRS4bin* drs4bin = new DRS4bin(ifname);
   drs4bin->Convert(entry1,entry2);
}
