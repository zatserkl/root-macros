//--            WET     b1 mean db1     sigma1  b0 mean db0     sigma0
//--------------------------------------------------------------------
//--DATA	0	1223	0.1	32.94	1029	0.1	30.53
//--DATA	13.32	1303	0.1	34.41	1061	0.1	34.49
//--DATA	26.54	1339	0.1	34.49	1150	0.1	41.32
//--DATA	39.87	1383	0.1	34.95	1337	0.2	58.07
//--DATA	52.91	1376	0.1	34.63	1556	0.1	38.87
//--DATA	79.45	1495	0.1	37.0	1313	0.1	38.51
//--DATA	105.82	1647	0.1	41.48	1020	0.1	41.28
//--DATA	132.36	1911	0.2	56.41	632	0.2	52.44
//--DATA	145.68	2156	0.4	76.43	344	0.3	69.9
//--DATA	150.72	2294	0.5	95.97	204.8	0.9	84.12
//--DATA	158.73	2557	0.2	38.43	87.23	0.63	18.08
//--DATA	185.27	2170	0.2	52.15	0	100	100
//--DATA	198.59	1938	0.2	52.41	0	100	100
//--DATA	211.63	1695	0.3	58.16	0	100	100
//--DATA	224.96	1401	0.2	70.41	0	100	100
//--DATA	238.18	1070	0.4	79.21	0	100	100
//--DATA	251.50	648.9	0.4	100.9	0	100	100
//--DATA	256.53	555.5	0.5	117.5	0	100	100
//--DATA	261.57	270	0.9	117.3	0	100	100
//--DATA	264.72	198.1	7.7	108.5	0	100	100
//--DATA	269.75	0	100	100	0	100	100

#include <TGraphErrors.h>
#include <TCanvas.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

using std::cout;     using std::endl;

void read_data()
{
   Float_t wet[100];
   Float_t b1[100];
   Float_t db1[100];
   Float_t s1[100];
   Float_t b0[100];
   Float_t db0[100];
   Float_t s0[100];
   Int_t np = 0;
   
   ifstream thisfile(__FILE__);
   if (!thisfile) {
      cout<< "File not found: " << __FILE__ <<endl;
      return;
   }
   std::string line;
   while (std::getline(thisfile,line))
   {
      const std::string pattern = "//--DATA";
      if (line.find(pattern) == 0)              // pattern position == 0
      {
         std::stringstream ss(line);
         ss.ignore(pattern.size());             // skip pattern
         ss >> wet[np] >> b1[np] >> db1[np] >> s1[np] >> b0[np] >> db0[np] >> s0[np];
         cout<< np << "\t wet " << wet[np] << " b1 " << b1[np] << " db1 " << db1[np] << " s1 " << s1[np] <<endl;
         np++;
      }
   }
   thisfile.close();

   // plot b1

   TGraphErrors* gr_b1 = new TGraphErrors(np, wet, b1, 0, db1);
   gr_b1->SetNameTitle("gr_b1", "Upstream bulky;WET, mm");
   gr_b1->SetMarkerStyle(20);
   gr_b1->SetMarkerColor(2);
   gr_b1->SetLineColor(2);

   new TCanvas;
   gr_b1->Draw("ap");
}
