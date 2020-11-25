#ifndef histogram
#define histogram

#include "histogram.C"

void histogram::set_bin_number(const int& n_bins){m_n_bins = n_bins; }
void histogram::set_variable(const std::string& variable){m_variable = variable; }
void histogram::set_x_limits(const float& x_min, const float& x_max){m_x_min = x_min; m_x_max = x_max; }
void histogram::set_y_limits(const float& y_min, const float& y_max){m_y_min = y_min; m_y_max = y_max; }
void histogram::set_x_label(const std::string& x_label){m_x_label = x_label; }
void histogram::set_y_label(const std::string& y_label){m_y_label = y_label; }
void histogram::set_weight(const string& weight){m_weight = weight; }
void histogram::set_stats(const bool& show_stats){
  if(show_stats == true){m_show_stats = 1; }
  else{m_show_stats = 1; }
}
void histogram::set_custom_x_range(const string& lower_axis_limit, const string& upper_axis_limit){
  m_plot_custom_axis = true;
  m_lower_axis_bounds = lower_axis_limit;
  m_upper_axis_bounds = upper_axis_limit;
}
void histogram::set_hist_title(const std::string& hist_title){m_hist_title = hist_title; }
void histogram::set_legend_names(const vector<std::string>& legend_names){ m_show_legend = true; m_legend_names = legend_names; }
void histogram::load_data_from_file(const std::string& list_of_files, const std::string& tree_name){
    std::vector<std::string> root_files = TxtFileListToVecStr(list_of_files);
    for(int i=0; i<root_files.size(); i++){data_files.push_back(ROOT::RDataFrame(tree_name, root_files[i])); std::cout << "Loaded: " << root_files[i] << std::endl; }
    set_hist_names(root_files);
}
void histogram::plot_hist(const std::string& filename){
  m_canvas = new TCanvas("c","c",10,10,700,900);
  TText T; T.SetTextFont(42); T.SetTextAlign(21);
  TString c_text = m_hist_title;

  std::vector<TH1D*> hists;
  for(int i{0}; i<data_files.size(); i++){
    ROOT::RDataFrame *temp_df = &data_files[i]; 
    if(m_weight != "None"){ 
      auto temp_hist = temp_df->Histo1D({m_hist_names[i].c_str(), "", m_n_bins, m_x_min, m_x_max}, m_variable, m_weight); 
      hists.push_back((TH1D*)temp_hist->Clone());  
      }
    else{ 
      auto temp_hist = temp_df->Histo1D({m_hist_names[i].c_str(), "", m_n_bins, m_x_min, m_x_max}, m_variable);
      hists.push_back((TH1D*)temp_hist->Clone());
    }      
  }

  hists[0]->GetXaxis()->SetTitle(m_x_label.c_str());
  hists[0]->GetYaxis()->SetTitle(m_y_label.c_str());
  hists[0]->SetMaximum(m_y_max);
  hists[0]->SetMinimum(m_y_min);
  hists[0]->SetStats(0);

  vector<int> colours = {600, 820, 400-6, 840, 920, 616, 860, 632, 432, 880, 416+2, 800, 900};
  
  for(int i=0; i<data_files.size(); i++){
    if(i < colours.size()){hists[i]->SetLineColor(colours[i]); }
    else{hists[i]->SetLineColor(colours[i-colours.size()]); std::cout << "That is a lot of lines you are plotting!" << std::endl;} 
  }

  hists[0]->SetLineWidth(2);
  hists[0]->Draw("HIST");
  T.DrawTextNDC(.5,.95,c_text);
  for(int i{1}; i<data_files.size(); i++){
      hists[i]->SetLineWidth(2);
      hists[i]->Draw("HIST SAME");
  }

  if(m_plot_custom_axis == true){
    TAxis* axis = hists[0]->GetXaxis();
    axis->SetNdivisions(-502);
    axis->ChangeLabel(1,-1,-1,-1,-1,-1, m_lower_axis_bounds);
    axis->ChangeLabel(-1,-1,-1,-1,-1,-1,m_upper_axis_bounds);
  }

  auto legend = new TLegend(0.9,0.9,0.5,0.7);
  if(m_show_legend == true){
    for(int i=0; i<data_files.size(); i++){
      legend->AddEntry(m_hist_names[i].c_str(),m_legend_names[i].c_str(),"l");
    }
    legend->Draw();
  }

  m_canvas->SetLeftMargin(0.15);
  m_canvas->Update();
  m_canvas->Modified();
  m_canvas->SaveAs(filename.c_str());
  std::cout << "Saved Histogram: " << filename << std::endl;
  delete m_canvas;
  delete legend;
}






#endif  
