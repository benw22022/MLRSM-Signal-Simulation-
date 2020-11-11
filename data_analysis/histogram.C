vector<string> TxtFileListToVecStr(const string f) {
  std::ifstream infile(f, std::ifstream::in);

  vector<string> v;
  string line("");
  while (getline(infile, line)) {

    line.erase(std::remove(line.begin(), line.end(), '\n'),  line.end());
    line.erase(std::remove(line.begin(), line.end(), '\r'),  line.end());
    line.erase(std::remove(line.begin(), line.end(), '\t'),  line.end());
      if(line.find_first_not_of(' ') == std::string::npos) continue; 
      if(line[line.find_first_not_of(' ')] == '#') continue; 
      v.push_back(line);
    }
  infile.close();

  return v;
}

class histogram{

  private:
    int m_n_bins;
    float m_x_min;
    float m_x_max;
    float m_y_min;
    float m_y_max;

    std::string m_variable;
    std::string m_weight;

    std::string m_hist_title;
    std::string m_x_label;
    std::string m_y_label;
    int m_show_stats;

    bool m_plot_custom_axis;
    std::string m_lower_axis_bounds;
    std::string m_upper_axis_bounds;
    
    std::vector<std::string> m_hist_names;
    std::string m_input_file;

    bool m_show_legend;
    std::vector<std::string> m_legend_names;

    bool m_set_log_scale;

    TCanvas *m_canvas;
    std::vector<ROOT::RDataFrame> m_data_files;

    void set_hist_names(const std::vector<std::string>& files){
        for(int i{0}; i<files.size(); i++){
            std::size_t found = files[i].find_last_of("/");
            std::string str = files[i].substr(found+1, files[i].size()-5);
            m_hist_names.push_back(str);
        }
    }

  public:
      //Constructors
      histogram() :
         m_n_bins{0},
         m_x_min{0},
         m_x_max{0},
         m_y_min{0},
         m_y_max{0},
         m_variable{"None"},
         m_weight{"None"},
         m_hist_title{"None"},
         m_x_label{"None"},
         m_y_label{"None"},
         m_show_stats{0},
         m_plot_custom_axis{false},
         m_lower_axis_bounds{"None"},
         m_upper_axis_bounds{"None"} ,
         m_hist_names{},
         m_input_file{"None"},
         m_show_legend{false},
         m_legend_names{""},
         m_set_log_scale{false}
        {std::cout << "default Constructor called" << std::endl;}

      histogram(int n_bins, std::string variable) :
         m_n_bins{n_bins},
         m_x_min{0},
         m_x_max{0},
         m_y_min{0},
         m_y_max{0},
         m_variable{variable},
         m_weight{"None"},
         m_hist_title{"None"},
         m_x_label{"None"},
         m_y_label{"None"},
         m_show_stats{0},
         m_plot_custom_axis{false},
         m_lower_axis_bounds{"None"},
         m_upper_axis_bounds{"None"},
         m_hist_names{},
         m_input_file{"None"},
         m_show_legend{false},
         m_legend_names{""},
         m_set_log_scale{false}
        {std::cout << "Bins Variable Constructor called" << std::endl;}

      histogram(int n_bins, std::string variable, std::string weight) :
         m_n_bins{n_bins},
         m_x_min{0},
         m_x_max{0},
         m_y_min{0},
         m_y_max{0},
         m_variable{variable},
         m_weight{weight},
         m_hist_title{"None"},
         m_x_label{"None"},
         m_y_label{"None"},
         m_show_stats{0},
         m_plot_custom_axis{false},
         m_lower_axis_bounds{"None"},
         m_upper_axis_bounds{"None"}, 
         m_hist_names{},
         m_input_file{"None"},
         m_show_legend{false},
         m_legend_names{""},
         m_set_log_scale{false}
        {std::cout << "Bins Variable Constructor called" << std::endl;}

      ~histogram() { std::cout << "destructor called" << std::endl; }
      // Member funcions
      void set_bin_number(const int&);
      void set_variable(const std::string&);
      void set_x_limits(const float&, const float&);
      void set_y_limits(const float&, const float&);
      void set_x_label(const std::string&);
      void set_y_label(const std::string&);
      void set_weight(const string&);
      void set_hist_title(const std::string&);
      void set_stats(const bool& );
      void set_custom_x_range(const string&, const string&);
      void set_legend_names(const std::vector<std::string>&);
      void set_legend_names(std::string&);
      void disable_legend();
      void disable_custom_axis();
      void set_log_scale(const bool&);
      void load_file_list(const std::string&);
      void load_data_from_file(const std::string& root_files, const std::string& tree_name);
      void plot_hist(const std::string&);
      std::string get_variable(){return m_variable;}
};

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
void histogram::disable_legend(){m_show_legend = false;}
void histogram::disable_custom_axis(){m_plot_custom_axis = false;}
void histogram::set_custom_x_range(const string& lower_axis_limit, const string& upper_axis_limit){ 
  m_plot_custom_axis = true;
  m_lower_axis_bounds = lower_axis_limit;
  m_upper_axis_bounds = upper_axis_limit;
}
void histogram::set_log_scale(const bool& set_log_scale){ m_set_log_scale = set_log_scale; }
void histogram::set_hist_title(const std::string& hist_title){m_hist_title = hist_title; }
void histogram::set_legend_names(const vector<std::string>& legend_names){ m_show_legend = true; m_legend_names = legend_names; }
void histogram::set_legend_names(std::string& legend_template){
  
  if(m_input_file != "None"){
    std::vector<std::string> root_files = TxtFileListToVecStr(m_input_file);
    m_show_legend = true;
    std::vector<std::string> flags;
    std::vector<std::vector<std::string>> flag_values;
    std::string tmp_str;
    bool begin_flag = false;

    // Find flags
    for(int i=0; i<legend_template.length(); i++){

      if(strncmp(&legend_template[i],"$",1) == 0 && begin_flag == false){
        begin_flag = true;
      }
      else if(strncmp(&legend_template[i],"$",1) == 0 && begin_flag == true){
        begin_flag = false;
        std::cout << "Found flag: " << tmp_str << std::endl;
        flags.push_back(tmp_str);
        tmp_str.clear();
      }
      else if(begin_flag == true){
        tmp_str = tmp_str + legend_template.at(i);
      }
    }  
    
    // Find the values of those flags
    tmp_str.clear();
    for(int i=0; i<root_files.size(); i++){
      std::vector<std::string> tmp_vec;
      std::cout << root_files[i] << std::endl;
      for(int j=0; j<flags.size(); j++){

        std::size_t found = root_files[i].find(flags[j]);

        std::cout << "found flag " << flags[j] << " at position " << found << std::endl;

        bool begin_flag = false;

        for(size_t k=found; k<root_files[i].length(); k++){
          
          std::string tmp_file = root_files[i];

          if(strncmp(&tmp_file[k],"_",1) == 0 && begin_flag == false){
            begin_flag = true;
          }
          else if(strncmp(&tmp_file[k],"_",1) == 0 && begin_flag == true){
            begin_flag = false;
            tmp_vec.push_back(tmp_str);
            tmp_str.clear();
            break;
          }
          else if(begin_flag == true){
            tmp_str = tmp_str + root_files[i].at(k);
          }
        }    
      }
      flag_values.push_back(tmp_vec);
    }

   m_legend_names.clear();
   for(int i=0; i<flag_values.size(); i++){ 
     std::string legend = "M_{W_{R}} = " + flag_values[i][0] + "      M_{N_{5}} = " + flag_values[i][1] + "     M_{N_{4,6}} = " + flag_values[i][2];
     std::cout << legend << std::endl;
     m_legend_names.push_back(legend);
   }
   
  }
  else{std::cout << "ERROR: Input file not set" << std::endl;}
}
void histogram::load_data_from_file(const std::string& list_of_files, const std::string& tree_name){
    m_input_file = list_of_files;
    std::vector<std::string> root_files = TxtFileListToVecStr(list_of_files);
    for(int i=0; i<root_files.size(); i++){m_data_files.push_back(ROOT::RDataFrame(tree_name, root_files[i])); std::cout << "Loaded: " << root_files[i] << std::endl; }
    set_hist_names(root_files);
}
void histogram::plot_hist(const std::string& filename){
  m_canvas = new TCanvas("c","c",10,10,700,900);
  TText T; T.SetTextFont(42); T.SetTextAlign(21);
  TString c_text = m_hist_title;

  std::vector<TH1D*> hists;
  for(int i{0}; i<m_data_files.size(); i++){
    ROOT::RDataFrame *temp_df = &m_data_files[i]; 
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

  vector<int> colours = {600, 840+5, 900, 920, 616, 860, 632,  432, 880, 416+2, 800,  820-2,  400-6};
  
  for(int i=0; i<m_data_files.size(); i++){
    if(i < colours.size()){hists[i]->SetLineColor(colours[i]); }
    else{hists[i]->SetLineColor(colours[i-colours.size()]); std::cout << "That is a lot of lines you are plotting!" << std::endl;} 
  }

  hists[0]->SetLineWidth(2);
  hists[0]->Draw("HIST");
  T.DrawTextNDC(.5,.95,c_text);
  for(int i{1}; i<m_data_files.size(); i++){
      hists[i]->SetLineWidth(2);
      hists[i]->Draw("HIST SAME");
  }

  if(m_set_log_scale == true){
    m_canvas->SetLogy();
  }


  TAxis* axis = hists[0]->GetXaxis();
  if(m_plot_custom_axis == true){
    axis->SetNdivisions(-502);
    axis->ChangeLabel(1,-1,-1,-1,-1,-1, m_lower_axis_bounds);
    axis->ChangeLabel(-1,-1,-1,-1,-1,-1,m_upper_axis_bounds);
  }

  auto legend = new TLegend(0.9,0.9,0.5,0.7);
  if(m_show_legend == true){
    for(int i=0; i<m_data_files.size(); i++){
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





