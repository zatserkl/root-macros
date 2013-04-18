void plot_endcapOccupancyMap(int run)
{
   TH2* endcapOccupancyMap = 0;

   std::stringstream ss_file;
   ss_file.str("");
   ss_file << "DQM_Pixel_R0000" << run << ".root";
   TFile* file = TFile::Open(ss_file.str().c_str());
   if (!file) {
      cout<< "File not found: " << file.str().c_str() <<endl;
      return;
   }

   std::stringstream ss_dir;
   ss_dir.str("");

   ss_dir << "/DQMData/Run " << run << "/Pixel/Run summary/Endcap";
   file->cd(ss_dir.str().c_str());

   endcapOccupancyMap = (TH2*) gDirectory->Get("endcapOccupancyMap");
   
   // new TCanvas;
   new TCanvas("can_endcapOccupancyMap", Form("endcapOccupancyMap for run %d",run), 800, 600);
   endcapOccupancyMap->Draw("lego2");

   //---------------------------------------------------------------
   
   // std::stringstream ss_file_dir;
   // ss_file_dir.str("");

   // ss_file_dir << "DQM_Pixel_R0000" << run << ".root:" << "/DQMData/Run " << run << "/Pixel/Run summary/Endcap";
   // gROOT->cd(ss_file_dir.str().c_str());

   // endcapOccupancyMap = (TH2*) gDirectory->Get("endcapOccupancyMap");
   // 
   // new TCanvas;
   // endcapOccupancyMap->Draw("lego2");
}
