#include <sys/types.h>
#include <sys/stat.h>
#include <ROOT/RDataFrame.hxx>
#include <TCanvas.h>
#include <vector>
#include <TLorentzVector.h>
#include <TFile.h>
#include <utility>     
#include <math.h>
#include "histogram.C"

int plot_distributions(){

    // Enable multithreading for RDataFrame
    ROOT::EnableImplicitMT();
	
    // Use sum of weights squared, makes sure error calculated properly
    TH1::SetDefaultSumw2();

    std::string file_extn = ".png";

    //----------------------------------------------------UNWEIGHTED----------------------------------------------------------------//
    float y_min = 0;
    float y_max = 1.0e-3;
    histogram hist(40, "eTmiss");
    hist.set_x_limits(0, 500);
    hist.set_y_limits(y_min, 3000);
    hist.set_x_label("E_{T}^{misss} (GeV)");
    hist.set_y_label("Events/Bin");
    hist.set_hist_title("eTmiss WRWR");
    hist.load_data_from_file("list_of_root_files.txt", "output_tree");
    //std::vector<std::string> legend_names = {"M_{W_{R}} = 1 TeV  \nM_{ N_{ 4}} = M_{ N_{ 5}} = M_{ N_{ 6}} = 2 TeV", "M_{W_{R}} = 1 TeV    \nM_{ N_{ 5}} = 2 TeV M_{ N_{ 4}} = M_{ N_{ 6}} = 10 TeV"};
    std::string legend_template = "M_{W} = $WR$ TeV    M_{N_{5} = $MN5$ TeV    M_{N_{4, 6}} = $MN46$ TeV";
    hist.set_legend_names(legend_template);
    //hist.set_log_scale(true);
    hist.plot_hist("../plots/unweighted/eTmiss_WRWR_plot"+file_extn);

    hist.set_variable("mjj");
    hist.set_x_limits(0, 3);
    hist.set_y_limits(y_min, 5000);
    hist.set_x_label("m_{jj} (TeV)");
    hist.set_hist_title("m_jj WRWR");
    hist.plot_hist("../plots/unweighted/mjj_WRWR_plot"+file_extn);

    hist.set_variable("jet1_pT");
    hist.set_x_limits(0, 2.5);
    hist.set_y_limits(y_min, 4000);
    hist.set_x_label("p_{T}^{j_{1}} (TeV)");
    hist.set_hist_title("jet1_pT WRWR");
    hist.plot_hist("../plots/unweighted/jet1_PT_WRWR_plot"+file_extn);

    hist.set_variable("jet2_pT");
    hist.set_x_limits(0, 1.4);
    hist.set_y_limits(y_min, 3000);
    hist.set_x_label("p_{T}^{j_{2}} (TeV)");
    hist.set_hist_title("jet2_pT WRWR");
    hist.plot_hist("../plots/unweighted/jet2_PT_WRWR_plot"+file_extn);
    
    hist.set_variable("delta_rap_jj");
    hist.set_x_limits(-5, 5);
    hist.set_y_limits(y_min, 1500);
    hist.set_x_label("#Deltay_{jj}");
    hist.set_hist_title("Delta_rap_jj WRWR");
    hist.plot_hist("../plots/unweighted/delta_rap_jj_WRWR_plot"+file_extn);

    hist.set_variable("dilep_centrality");
    hist.set_x_limits(0, 15);
    hist.set_y_limits(y_min, 4000);
    hist.set_x_label("|z_{ll}^{*}|");
    hist.set_hist_title("dilep_centrality WRWR");
    hist.plot_hist("../plots/unweighted/dilep_centrality_WRWR_plot"+file_extn);

    //----------------------------------------------------WEIGHTED----------------------------------------------------------------//
    hist.set_weight("weight");
    hist.set_variable("eTmiss");
    hist.set_x_limits(0, 500);
    hist.set_y_limits(0, 0.1);
    hist.set_x_label("E_{T}^{misss} (GeV)");
    hist.set_y_label("Events/Bin");
    hist.set_hist_title("eTmiss WRWR");
    hist.plot_hist("../plots/weighted/eTmiss_WRWR_plot"+file_extn);

    hist.set_variable("mjj");
    hist.set_x_limits(0, 3);
    hist.set_y_limits(0, 0.1);
    hist.set_x_label("m_{jj} (TeV)");
    hist.set_hist_title("m_jj WRWR");
    hist.plot_hist("../plots/weighted/mjj_WRWR_plot"+file_extn);

    hist.set_variable("jet1_pT");
    hist.set_x_limits(0, 2.5);
    hist.set_y_limits(0, 0.1);
    hist.set_x_label("p_{T}^{j_{1}} (TeV)");
    hist.set_hist_title("jet1_pT WRWR");
    hist.plot_hist("../plots/weighted/jet1_PT_WRWR_plot"+file_extn);

    hist.set_variable("jet2_pT");
    hist.set_x_limits(0, 1.4);
    hist.set_y_limits(0, 0.1);
    hist.set_x_label("p_{T}^{j_{2}} (TeV)");
    hist.set_hist_title("jet2_pT WRWR");
    hist.plot_hist("../plots/weighted/jet2_PT_WRWR_plot"+file_extn);
    
    hist.set_variable("delta_rap_jj");
    hist.set_x_limits(-5, 5);
    hist.set_y_limits(0, 0.1);
    hist.set_x_label("#Deltay_{jj}");
    hist.set_hist_title("Delta_rap_jj WRWR");
    hist.plot_hist("../plots/weighted/delta_rap_jj_WRWR_plot"+file_extn);

    hist.set_variable("dilep_centrality");
    hist.set_x_limits(0, 15);
    hist.set_y_limits(0, 0.1);
    hist.set_x_label("|z_{ll}^{*}|");
    hist.set_hist_title("dilep_centrality WRWR");
    hist.plot_hist("../plots/weighted/dilep_centrality_WRWR_plot"+file_extn);

    hist.set_variable("lep1_centrality");
    hist.set_x_limits(0, 15);
    hist.set_y_limits(0, 0.1);
    hist.set_x_label("|z_{l1}^{*}|");
    hist.set_hist_title("lep1_centrality WRWR");
    hist.plot_hist("../plots/weighted/lep1_centrality_WRWR_plot"+file_extn);

return 0;
}