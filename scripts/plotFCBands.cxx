void SetLHCbStyle(std::string);
#include "tools.h"
#include <stdlib.h> 
// ===========================================================================
std::pair<double,double> GetLimits(double,double,double,double,double,double,double,double);
double Getx1(double,double,double);
double Getx2(double,double,double);

void plotFCBands(std::string mode = "DsPhi"){


    std::string name;
    if (mode == "DsPhi") name = "#it{B^{+}#rightarrowD_{s}^{+}#phi}";
    if (mode == "DKst0") name = "it{B^{+}#rightarrowD^{+}K^{*0}}";

    std::cout << "Using mode name: " << mode << " Name: " << name << std::endl;
    SetLHCbStyle("norm");

    std::string sensitivityDir = "sensitivityDir/";
    std::vector<std::string> fileName;
    getdir(sensitivityDir,fileName);
    std::map<double,int> gen_count;
    std::map<double,std::vector<double>> map_gen_br_fit_br;
    std::map<double,std::vector<double>> map_gen_br_fit_stat;
    std::map<double,std::vector<double>> map_gen_br_fit_stat_hi;
    std::map<double,std::vector<double>> map_gen_br_fit_stat_lo;

    int n_points = 0;
    for(std::vector<std::string>::iterator it=fileName.begin();it!=fileName.end();it++){
        if( (*it).find("sensitivity")  ==std::string::npos ) continue;
        n_points++;
    }

    double* n_gen_br       = new double[n_points];
    double* n_fit_br       = new double[n_points];
    double* n_fit_stat     = new double[n_points];
    double* n_fit_stat_hi  = new double[n_points];
    double* n_fit_stat_lo  = new double[n_points];

    int n_points_new = 0;
    int count = 0;
    int nFits=0;
    int nBad=0;
    int nFPD = 0;
    int nBadMINOS = 0;
    for(std::vector<std::string>::iterator it=fileName.begin();it!=fileName.end();it++){
        if( (*it).find("sensitivity")  ==std::string::npos ) continue;
        if(count%1000==0) std::cout << "Filename: " << *it << std::endl;
    
        std::vector<std::string> tokens;
        split((*it),tokens,"_");

        double br = atof(tokens[1].c_str());
        int seed  = atoi(tokens[2].c_str());
        int num   = atoi(tokens[3].c_str());


        TFile newfile((sensitivityDir+*it).c_str(),"READ");
        if( !(newfile.GetListOfKeys()->Contains("toy_withSig"))  ) continue; 
        RooFitResult* result = (RooFitResult*)newfile.Get("toy_withSig");
        if(result->covQual()<2){continue;}
        if(result->covQual()<3){continue;}
 
        RooRealVar* param_final = (RooRealVar*) (result->floatParsFinal()).find("Branching_fraction");
        double fitted_br   = param_final->getVal();
        double fitted_stat = param_final->getError();
        double fitted_stat_hi = param_final->getAsymErrorHi();
        double fitted_stat_lo = param_final->getAsymErrorLo();
        if(count%1000==0)  std::cout << "Br: " << br << " Fitted Br: " << fitted_br << " +/- " << fitted_stat <<" (+" <<  fitted_stat_hi <<"," <<  fitted_stat_lo << ") seed: " << seed << " Number: " << num << std::endl;
        n_gen_br[count]      = br;
        n_fit_br[count]      = fitted_br;
        n_fit_stat[count]    = fitted_stat;
        n_fit_stat_hi[count] = fitted_stat_hi;
        n_fit_stat_lo[count] = fitted_stat_lo;
        count++ ;
        n_points_new++;
        gen_count[br]++;
        map_gen_br_fit_br[br].push_back(fitted_br);
        map_gen_br_fit_stat[br].push_back(fitted_stat);
        map_gen_br_fit_stat_hi[br].push_back(fitted_stat_hi);
        map_gen_br_fit_stat_lo[br].push_back(fitted_stat_lo);

    }

    int n_gen_points = map_gen_br_fit_br.size();
    
    double* gen_br        = new double[n_gen_points];
    double* fit_br        = new double[n_gen_points];
    double* fit_stat      = new double[n_gen_points];
    double* fit_stat_hi   = new double[n_gen_points];
    double* fit_stat_lo   = new double[n_gen_points];


   /*
    std::cout << "Gen points: " << n_gen_points <<std::endl;
    double* n_gen_line_br   = new double[n_gen_points];
    double* n_fit_line_br   = new double[n_gen_points];
    
    double* n_gen_line2_br   = new double[n_gen_points];
    double* n_fit_line2_br   = new double[n_gen_points];
    
    double* n_gen_line_br_lower   = new double[n_gen_points];
    double* n_fit_line_br_lower   = new double[n_gen_points];
    
    double* n_gen_line2_br_lower   = new double[n_gen_points];
    double* n_fit_line2_br_lower   = new double[n_gen_points];
    
    double* n_gen_line_br_sys   = new double[n_gen_points];
    double* n_fit_line_br_sys   = new double[n_gen_points];
    
    double* n_gen_line2_br_sys   = new double[n_gen_points];
    double* n_fit_line2_br_sys   = new double[n_gen_points];
    
    double* n_gen_line_br_sys_lower   = new double[n_gen_points];
    double* n_fit_line_br_sys_lower   = new double[n_gen_points];
    
    double* n_gen_line2_br_sys_lower   = new double[n_gen_points];
    double* n_fit_line2_br_sys_lower   = new double[n_gen_points];


    // FC arrays 
    double* FC_lower_line_95   = new double[n_gen_points];
    double* FC_upper_line_95   = new double[n_gen_points];


    
    
    */

    // First plot Nfit vs. Ngen and Sigmafit vs. Ngen 
    count = 0;
    for(std::map<double,std::vector<double>>::iterator it=map_gen_br_fit_br.begin();it!=map_gen_br_fit_br.end();it++){ 
        std::vector<double> temp = it->second;
        double gen_value = it->first;
        std::sort(temp.begin(), temp.end());
        std::vector<double> temp_stat    =  map_gen_br_fit_stat[gen_value];
        std::vector<double> temp_stat_hi =  map_gen_br_fit_stat_hi[gen_value];
        std::vector<double> temp_stat_lo =  map_gen_br_fit_stat_lo[gen_value];

        //TCanvas *cav = new TCanvas(Form("Canvas_%d",count),Form("Distribution_%d",count),200,10,1200,500);
        //cav->Divide(4,1);
        //cav->cd(1);
        
        // Draw fitted value
        TH1D *h_temp = new TH1D(Form("h_fit_%d",count), "", 100,-10.0,20.0);
        for(int i = 0; i < temp.size(); i++){
          h_temp->Fill(temp[i]);
        }
        h_temp->Draw();
        
        TF1 *gaussian = new TF1("Gaussian","gaus",-10.0,20);
        h_temp->Fit("Gaussian","Q0R");
        double mean  = gaussian->GetParameter(1);
        double sigma = gaussian->GetParameter(2);
        const double* err=gaussian->GetParErrors();
        double meanerr=err[1];
        double sigmaerr=err[2];
        gaussian->SetLineColor(kBlue); 
        gaussian->Draw("same");
        TPaveLabel *pav1 = new TPaveLabel(0.5,0.89,0.93,0.94,Form("#bf{Mean: } %f #pm %f", mean,  meanerr),"NDC");
        TPaveLabel *pav2 = new TPaveLabel(0.5,0.81,0.93,0.86,Form("#bf{Sigma:} %f #pm %f",sigma, sigmaerr),"NDC");
        pav1->SetTextColor(kRed); pav2->SetTextColor(kRed);
        pav1->SetBorderSize(0);   pav2->SetBorderSize(0);
        //pav1->SetFillStyle(4000); pav2->SetFillStyle(4000);
        pav1->SetFillColor(0);    pav2->SetFillColor(0);
        pav1->Draw();             pav2->Draw();

        // Draw stat error
        //cav->cd(2);
        TH1D *h_temp_stat = new TH1D(Form("h_fit_stat_%d",count), "", 100,0.0,5.0);
        for(int i = 0; i < temp_stat.size(); i++){
          h_temp_stat->Fill(temp_stat[i]);
        }
        h_temp_stat->Draw();  

        TF1 *gaussian2 = new TF1("Gaussian","gaus",-10.0,20);
        h_temp_stat->Fit("Gaussian","Q0R");
        double mean_stat  = gaussian2->GetParameter(1);
        double sigma_stat = gaussian2->GetParameter(2);
        err=gaussian2->GetParErrors();
        meanerr=err[1];
        sigmaerr=err[2];
        gaussian2->SetLineColor(kBlue); 
        gaussian2->Draw("same");
        TPaveLabel *pav3 = new TPaveLabel(0.5,0.89,0.93,0.94,Form("#bf{Mean: } %f #pm %f", mean_stat,  meanerr),"NDC");
        TPaveLabel *pav4 = new TPaveLabel(0.5,0.81,0.93,0.86,Form("#bf{Sigma:} %f #pm %f",sigma_stat, sigmaerr),"NDC");
        pav3->SetTextColor(kRed); pav4->SetTextColor(kRed);
        pav3->SetBorderSize(0);   pav4->SetBorderSize(0);
        //pav3->SetFillStyle(4000); pav4->SetFillStyle(4000);
        pav3->SetFillColor(0);    pav4->SetFillColor(0);
        pav3->Draw();             pav4->Draw();      
        
        // Draw stat error hi
        //cav->cd(3);
        TH1D *h_temp_stat_hi = new TH1D(Form("h_fit_stat_hi_%d",count), "", 100,-1.0,5.0);
        for(int i = 0; i < temp_stat_hi.size(); i++){
          h_temp_stat_hi->Fill(temp_stat_hi[i]);
        }
        h_temp_stat_hi->Draw();     

        TF1 *gaussian3 = new TF1("Gaussian","gaus",-1.0,5.0);
        h_temp_stat_hi->Fit("Gaussian","Q0R");
        double mean_stat_hi  = gaussian3->GetParameter(1);
        double sigma_stat_hi = gaussian3->GetParameter(2);
        err=gaussian3->GetParErrors();
        meanerr=err[1];
        sigmaerr=err[2];
        gaussian3->SetLineColor(kBlue); 
        gaussian3->Draw("same");
        TPaveLabel *pav5 = new TPaveLabel(0.5,0.89,0.93,0.94,Form("#bf{Mean: } %f #pm %f", mean_stat_hi,  meanerr),"NDC");
        TPaveLabel *pav6 = new TPaveLabel(0.5,0.81,0.93,0.86,Form("#bf{Sigma:} %f #pm %f",sigma_stat_hi, sigmaerr),"NDC");
        pav5->SetTextColor(kRed); pav6->SetTextColor(kRed);
        pav5->SetBorderSize(0);   pav6->SetBorderSize(0);
        //pav5->SetFillStyle(4000); pav6->SetFillStyle(4000);
        pav5->SetFillColor(0);    pav6->SetFillColor(0);
        pav5->Draw();             pav6->Draw();    
        
        // Draw stat error lo
        //cav->cd(4);
        TH1D *h_temp_stat_lo = new TH1D(Form("h_fit_stat_lo_%d",count), "", 100,-1.0,5.0);
        for(int i = 0; i < temp_stat_lo.size(); i++){
          h_temp_stat_lo->Fill(-temp_stat_lo[i]);
        }
        h_temp_stat_lo->Draw();

        TF1 *gaussian4 = new TF1("Gaussian","gaus",-1.0,5.0);
        h_temp_stat_lo->Fit("Gaussian","Q0R");
        double mean_stat_lo  = gaussian4->GetParameter(1);
        double sigma_stat_lo = gaussian4->GetParameter(2);
        err=gaussian4->GetParErrors();
        meanerr=err[1];
        sigmaerr=err[2];
        gaussian4->SetLineColor(kBlue); 
        gaussian4->Draw("same");
        TPaveLabel *pav7 = new TPaveLabel(0.5,0.89,0.93,0.94,Form("#bf{Mean: } %f #pm %f", mean_stat_lo,  meanerr),"NDC");
        TPaveLabel *pav8 = new TPaveLabel(0.5,0.81,0.93,0.86,Form("#bf{Sigma:} %f #pm %f",sigma_stat_lo, sigmaerr),"NDC");
        pav7->SetTextColor(kRed); pav8->SetTextColor(kRed);
        pav7->SetBorderSize(0);   pav8->SetBorderSize(0);
        pav7->SetFillColor(0);    pav8->SetFillColor(0);

        //pav7->SetFillStyle(4000); pav8->SetFillStyle(4000);
        pav7->Draw();             pav8->Draw();  

        std::cout << "Gen BR: " << gen_value << " Length: " << temp.size() << "\t Mean: " << mean<< "\tStat: " << mean_stat<<"\t(+" << mean_stat_hi<< ",-"<< mean_stat_lo<<")" << std::endl;

        //Add means to arrays
        gen_br[count]        = gen_value;
        fit_br[count]        = mean;
        fit_stat[count]      = mean_stat*mean_stat;
        fit_stat_hi[count]   = mean_stat_hi*mean_stat_hi;
        fit_stat_lo[count]   = mean_stat_lo*mean_stat_lo;

        count++;
    }
        
    //TCanvas *cav2 = new TCanvas("Canvas_trends","Distribution_trends",200,10,1200,500);
    //cav2->Divide(4,1);
    
    //cav2->cd(1);
    TGraph* gr_fit_br = new TGraph(n_gen_points, gen_br, fit_br);
    gr_fit_br->Draw("AP");
    TF1 *fitline_br = new TF1("fitline_br","[0]+[1]*x",0.0,8.0); 
    gr_fit_br->Fit("fitline_br","Q0R");
    double p0_mean = fitline_br->GetParameter(0);
    double p1_mean = fitline_br->GetParameter(1);
    fitline_br->SetLineColor(kRed);
    fitline_br->Draw("same");
    err=fitline_br->GetParErrors();
    double p0_meanerr=err[0];
    double p1_meanerr=err[1];
    std::cout << "p0_mean:  " << p0_mean << "\t+/- " << p0_meanerr <<std::endl;
    std::cout << "p1_mean:  " << p1_mean << "\t+/- " << p1_meanerr <<std::endl;

    //cav2->cd(2);
    TGraph* gr_fit_stat = new TGraph(n_gen_points, gen_br, fit_stat);
    gr_fit_stat->Draw("AP");
    TF1 *fitline_stat = new TF1("fitline_stat","[0]+[1]*x",0.0,8.0); 
    gr_fit_stat->Fit("fitline_stat","Q0R");
    double p0_stat = fitline_stat->GetParameter(0);
    double p1_stat = fitline_stat->GetParameter(1);
    fitline_stat->SetLineColor(kRed);
    fitline_stat->Draw("same");
    err=fitline_stat->GetParErrors();
    double p0_staterr=err[0];
    double p1_staterr=err[1];
    std::cout << "p0_stat:  " << p0_stat << "\t+/- " << p0_staterr <<std::endl;
    std::cout << "p1_stat:  " << p1_stat << "\t+/- " << p1_staterr <<std::endl;


    //cav2->cd(3);
    TGraph* gr_fit_stat_hi = new TGraph(n_gen_points, gen_br, fit_stat_hi);
    gr_fit_stat_hi->Draw("AP");    
    TF1 *fitline_stat_hi = new TF1("fitline_stat_hi","[0]+[1]*x",0.0,8.0); 
    gr_fit_stat_hi->Fit("fitline_stat_hi","Q0R");
    double p0_stat_hi = fitline_stat_hi->GetParameter(0);
    double p1_stat_hi = fitline_stat_hi->GetParameter(1);
    fitline_stat_hi->SetLineColor(kRed);
    fitline_stat_hi->Draw("same");
    err=fitline_stat_hi->GetParErrors();
    double p0_stat_hierr=err[0];
    double p1_stat_hierr=err[1];
    std::cout << "p0_stat_hi:  " << p0_stat_hi << "\t+/- " << p0_stat_hierr <<std::endl;
    std::cout << "p1_stat_hi:  " << p1_stat_hi << "\t+/- " << p1_stat_hierr <<std::endl;


    //cav2->cd(4);
    TGraph* gr_fit_stat_lo = new TGraph(n_gen_points, gen_br, fit_stat_lo);
    gr_fit_stat_lo->Draw("AP");
    TF1 *fitline_stat_lo = new TF1("fitline_stat_lo","[0]+[1]*x",0.0,8.0); 
    gr_fit_stat_lo->Fit("fitline_stat_lo","Q0R");
    double p0_stat_lo = fitline_stat_lo->GetParameter(0);
    double p1_stat_lo = fitline_stat_lo->GetParameter(1);
    fitline_stat_lo->SetLineColor(kRed);
    fitline_stat_lo->Draw("same");
    err=fitline_stat_lo->GetParErrors();
    double p0_stat_loerr=err[0];
    double p1_stat_loerr=err[1];
    std::cout << "p0_stat_lo:  " << p0_stat_lo << "\t+/- " << p0_stat_loerr <<std::endl;
    std::cout << "p1_stat_lo:  " << p1_stat_lo << "\t+/- " << p1_stat_loerr <<std::endl;

    double p0_sys = 0.78;
    double p1_sys = 0.0;
    double CL = 0.95;

    //std::pair<double,double> limits1 = GetLimits(2.0,p0_mean,p1_mean,p0_stat,p1_stat,p0_sys,p1_sys,CL);
    //std::pair<double,double> limits2 = GetLimits(4.0,p0_mean,p1_mean,p0_stat,p1_stat,p0_sys,p1_sys,CL);
    //std::pair<double,double> limits3 = GetLimits(6.0,p0_mean,p1_mean,p0_stat,p1_stat,p0_sys,p1_sys,CL);

    int n_plot_points = 20;
    double gen_min = 0.1;
    double gen_max = 7.0;

    double* FC_gen           = new double[n_plot_points];
    double* FC_95_stat_lower = new double[n_plot_points];
    double* FC_95_stat_upper = new double[n_plot_points];
    double* FC_95_syst_lower = new double[n_plot_points];
    double* FC_95_syst_upper = new double[n_plot_points];
    double* FC_90_stat_lower = new double[n_plot_points];
    double* FC_90_stat_upper = new double[n_plot_points];
    double* FC_90_syst_lower = new double[n_plot_points];
    double* FC_90_syst_upper = new double[n_plot_points];
    count = 0;

    for(int i = 0; i<n_plot_points; i++){
        double gen_br_temp = gen_min + i*(gen_max-gen_min)/n_plot_points;
        FC_gen[count] = gen_br_temp;
        std::pair<double,double> limits_90_stat = GetLimits(gen_br_temp,p0_mean,p1_mean,p0_stat,p1_stat,0.0,   0.0,   0.90);
        std::pair<double,double> limits_95_stat = GetLimits(gen_br_temp,p0_mean,p1_mean,p0_stat,p1_stat,0.0,   0.0,   0.95);
        std::pair<double,double> limits_95_syst = GetLimits(gen_br_temp,p0_mean,p1_mean,p0_stat,p1_stat,p0_sys,p1_sys,0.95);

        FC_95_stat_lower[count] = limits_95_stat.first;
        FC_95_stat_upper[count] = limits_95_stat.second;
        
        FC_90_stat_lower[count] = limits_90_stat.first;
        FC_90_stat_upper[count] = limits_90_stat.second;
        
        FC_95_syst_lower[count] = limits_95_syst.first;
        FC_95_syst_upper[count] = limits_95_syst.second;
        
        std::cout << "Gen BR: " << gen_br_temp << "\t(" <<  limits_95_stat.first << "," << limits_95_stat.second<<")" <<std::endl;

        count++; 
    }
    
    TCanvas *cav = new TCanvas("Canvas_FC","Canvas_FC",200,10,700,500);


    TGraph* gr_95_stat_lower  = new TGraph(n_plot_points, FC_95_stat_lower, FC_gen);
    TGraph* gr_95_stat_upper  = new TGraph(n_plot_points, FC_95_stat_upper, FC_gen);

    TGraph* gr_90_stat_lower  = new TGraph(n_plot_points, FC_90_stat_lower, FC_gen);
    TGraph* gr_90_stat_upper  = new TGraph(n_plot_points, FC_90_stat_upper, FC_gen);

    TGraph* gr_95_syst_lower  = new TGraph(n_plot_points, FC_95_syst_lower, FC_gen);
    TGraph* gr_95_syst_upper  = new TGraph(n_plot_points, FC_95_syst_upper, FC_gen);

    gr_95_syst_lower->SetLineStyle(kDashed);
    gr_95_syst_upper->SetLineStyle(kDashed);


    TGraph *grshade_95 = new TGraph(2*n_plot_points);
    for (int i=0;i<n_plot_points;i++) {
      grshade_95->SetPoint(i,FC_95_stat_upper[i],FC_gen[i]);
      grshade_95->SetPoint(n_plot_points+i,FC_95_stat_lower[n_plot_points-i-1],FC_gen[n_plot_points-i-1]);
    }
    grshade_95->SetFillColor(kYellow);

    TGraph *grshade_90 = new TGraph(2*n_plot_points);
    for (int i=0;i<n_plot_points;i++) {
      grshade_90->SetPoint(i,FC_90_stat_upper[i],FC_gen[i]);
      grshade_90->SetPoint(n_plot_points+i,FC_90_stat_lower[n_plot_points-i-1],FC_gen[n_plot_points-i-1]);
    }
    grshade_90->SetFillColor(kGreen);

    TMultiGraph *mg = new TMultiGraph();
    mg->SetTitle(Form("Sensitivity %s;Fitted #it{B}(%s) (#times10^{-7});Generated #it{B}(%s) (#times10^{-7})",name.c_str(),name.c_str(),name.c_str()));

    mg->Add(grshade_95, "F");
    mg->Add(grshade_90, "F");
    mg->Add(gr_90_stat_lower, "L");
    mg->Add(gr_90_stat_upper, "L");
    mg->Add(gr_95_stat_lower, "L");
    mg->Add(gr_95_stat_upper, "L");
    mg->Add(gr_95_syst_lower, "L");
    mg->Add(gr_95_syst_upper, "L");
    mg->Draw("A");

    mg->GetXaxis()->SetRangeUser(-10.0,15.0);
    mg->GetYaxis()->SetRangeUser(0.5,7.0);


    TLine *line_diagonal = new TLine(0.0,0.0,7.0,7.0); 
    line_diagonal->SetLineStyle(kDashed);
    line_diagonal->SetLineWidth(2);
    line_diagonal->SetLineColor(kWhite);
    line_diagonal->Draw();

    //double limit     = 4.45;
    //double limit_sys = 4.85;

    //double limit     = 3.8;
    //double limit_sys = 4.2;

    double limit     = 4.6;
    double limit_sys = 4.9;

    TLine *line = new TLine(1.16,0.1,1.16,limit_sys);
    line->SetLineWidth(2);
    line->SetLineColor(kRed);
    line->Draw();    

    TLine *line2 = new TLine(-10.0,limit,1.16,limit);
    line2->SetLineStyle(kDashed);
    line2->SetLineWidth(2);
    line2->SetLineColor(kRed);
    line2->Draw();


    TLine *line4 = new TLine(-10.0,limit_sys,1.16,limit_sys);
    line4->SetLineStyle(kDashed);
    line4->SetLineWidth(2);
    line4->SetLineColor(kRed);
    line4->Draw();



    gPad->RedrawAxis();

    TLegend* leg = new TLegend(0.65,0.2,0.92,0.4);
    //TLegend* leg = new TLegend(0.2,0.44,0.5,0.65);
    leg->AddEntry(grshade_90,"90% CL Band","f");
    leg->AddEntry(grshade_95,"95% CL Band","f");
    leg->AddEntry(gr_95_syst_lower,"With Systematics","l");
    leg->Draw();

        /*  
        std::cout << "Gen BR: " << gen_value << " Length: " << temp.size()<< " Limit point: " << temp.size()*0.05 << " Rounded: " << floor(temp.size()*0.05 + 0.5) << " Value: " << temp[floor(temp.size()*0.05 + 0.5)-1] << std::endl;
        
        //std::cout << "Values: [" ; for(std::vector<double>::iterator itv=temp.begin();itv!=temp.end();itv++){std::cout << *itv <<", ";} std::cout << "]"<< std::endl;
        
        TCanvas *cav = new TCanvas(Form("Canvas_%d",count),Form("Distribution_%d",count),200,10,700,500);
        TH1D *h_temp = new TH1D(Form("h_%d",count), "", 100,-10.0,20.0);
        for(int i = 0; i < temp.size(); i++){
          h_temp->Fill(temp[i]);
        }
        h_temp->Draw();


        TF1 *gaussian = new TF1("Gaussian","gaus",-10.0,20);
        h_temp->Fit("Gaussian","Q0R");
        double mean  = gaussian->GetParameter(1);
        double sigma = gaussian->GetParameter(2);
        const double* err=gaussian->GetParErrors();
        double meanerr=err[1];
        double sigmaerr=err[2];
        gaussian->SetLineColor(kBlue); gaussian->Draw("same");

        double sigma_FC_sys = 0.788;
        double sigma_FC_total = sqrt(pow(sigma,2) + pow(sigma_FC_sys,2));
        
        double xmin = -60.0;
        double xmax =  20.0;
        double stepsize = 0.01;

        TF1 *f1 = new TF1("f1" , "(x>0)*(pow(x-[0],2)/pow([1],2)) + (x<0)*((pow([0],2) - 2*[0]*x)/pow([1],2))",xmin,xmax);

        f1->SetParameter(0,mean);
        f1->SetParameter(1,sigma);
        f1->Draw("same");

        TF1 *f2 = new TF1("f2" , "(x>0)*exp(-0.5*(pow(x-[0],2)/pow([1],2))) + (x<0)*exp(-0.5*((pow([0],2) - 2*[0]*x)/pow([1],2)))",xmin,xmax);

        f2->SetParameter(0,mean);
        f2->SetParameter(1,sigma);
        f2->Draw("same");
        
        TF1 *f1_sys = new TF1("f1_sys" , "(x>0)*(pow(x-[0],2)/pow([1],2)) + (x<0)*((pow([0],2) - 2*[0]*x)/pow([1],2))",xmin,xmax);

        f1_sys->SetParameter(0,mean);
        f1_sys->SetParameter(1,sigma_FC_total);
        f1_sys->Draw("same");

        TF1 *f2_sys = new TF1("f2_sys" , "(x>0)*exp(-0.5*(pow(x-[0],2)/pow([1],2))) + (x<0)*exp(-0.5*((pow([0],2) - 2*[0]*x)/pow([1],2)))",xmin,xmax);

        f2_sys->SetParameter(0,mean);
        f2_sys->SetParameter(1,sigma_FC_total);
        f2_sys->Draw("same");



        double final_x1_95 = xmin;
        double final_x2_95 = xmax;
        double integral_95 = 0.0;
        bool found_interval_95 = false;

        double final_x1_90 = xmin;
        double final_x2_90 = xmax;
        double integral_90 = 0.0;
        bool found_interval_90 = false;


        double final_x1_95_sys = xmin;
        double final_x2_95_sys = xmax;
        double integral_95_sys = 0.0;
        bool found_interval_95_sys = false;

        double final_x1_90_sys = xmin;
        double final_x2_90_sys = xmax;
        double integral_90_sys = 0.0;
        bool found_interval_90_sys = false;

        // Statistical only
        for(int i = 0; i< (floor((mean - xmin)/stepsize)); i++ ) {
            double x1  = i*stepsize + xmin;
            double Rx1 = f2->Eval(i*stepsize+ xmin);


            double x2  = xmax;
            double Rx2 = f2->Eval(xmax);
            bool found_x2   = false;

            for(int j = 0; j < floor(( xmax - mean)/stepsize); j++){
                //std::cout << "\tx2:    " << mean + j*stepsize << "\tR(x2): " << f2->Eval(mean + j*stepsize) << std::endl;
                if(f2->Eval(mean + j*stepsize)<Rx1 && !found_x2){
                    x2  = mean + j*stepsize;
                    Rx2 = f2->Eval(mean + j*stepsize);
                    found_x2 = true;
                }
            }

            //std::cout << "x1:    " <<  x1 << "\tR(x1): " << Rx1 << "\t x2:    " <<  x2 << "\tR(x2): " << Rx2 << "\tIntegral: "<< f2->Integral(x1,x2)/f2->Integral(xmin,xmax)<< std::endl;
            //std::cout << std::endl;

            if(f2->Integral(x1,x2)/f2->Integral(xmin,xmax)<0.95&&!found_interval_95){
                final_x1_95 = x1;
                final_x2_95 = x2;
                found_interval_95 = true;
                integral_95 = f2->Integral(x1,x2)/f2->Integral(xmin,xmax);
            }
            if(f2->Integral(x1,x2)/f2->Integral(xmin,xmax)<0.90&&!found_interval_90){
                final_x1_90 = x1;
                final_x2_90 = x2;
                found_interval_90 = true;
                integral_90 = f2->Integral(x1,x2)/f2->Integral(xmin,xmax);
            }

        }

        std::cout << "x1:    " <<  final_x1_95 << "\tR(x1): " << f2->Eval(final_x1_95) << std::endl;
        std::cout << "x2:    " <<  final_x2_95 << "\tR(x2): " << f2->Eval(final_x2_95) << std::endl;
        std::cout <<"Integral: "<< integral_95<< std::endl;

        std::cout << "x1:    " <<  final_x1_90 << "\tR(x1): " << f2->Eval(final_x1_90) << std::endl;
        std::cout << "x2:    " <<  final_x2_90 << "\tR(x2): " << f2->Eval(final_x2_90) << std::endl;
        std::cout <<"Integral: "<< integral_90<< std::endl;

        // With systematics
        for(int i = 0; i< (floor((mean - xmin)/stepsize)); i++ ) {
            double x1  = i*stepsize + xmin;
            double Rx1 = f2_sys->Eval(i*stepsize+ xmin);


            double x2  = xmax;
            double Rx2 = f2_sys->Eval(xmax);
            bool found_x2   = false;

            for(int j = 0; j < floor(( xmax - mean)/stepsize); j++){
                //std::cout << "\tx2:    " << mean + j*stepsize << "\tR(x2): " << f2_sys->Eval(mean + j*stepsize) << std::endl;
                if(f2_sys->Eval(mean + j*stepsize)<Rx1 && !found_x2){
                    x2  = mean + j*stepsize;
                    Rx2 = f2_sys->Eval(mean + j*stepsize);
                    found_x2 = true;
                }
            }

            //std::cout << "x1:    " <<  x1 << "\tR(x1): " << Rx1 << "\t x2:    " <<  x2 << "\tR(x2): " << Rx2 << "\tIntegral: "<< f2_sys->Integral(x1,x2)/f2_sys->Integral(xmin,xmax)<< std::endl;
            //std::cout << std::endl;

            if(f2_sys->Integral(x1,x2)/f2_sys->Integral(xmin,xmax)<0.95&&!found_interval_95_sys){
                final_x1_95_sys = x1;
                final_x2_95_sys = x2;
                found_interval_95_sys = true;
                integral_95_sys = f2_sys->Integral(x1,x2)/f2_sys->Integral(xmin,xmax);
            }
            if(f2_sys->Integral(x1,x2)/f2_sys->Integral(xmin,xmax)<0.90&&!found_interval_90_sys){
                final_x1_90_sys = x1;
                final_x2_90_sys = x2;
                found_interval_90_sys = true;
                integral_90_sys = f2_sys->Integral(x1,x2)/f2_sys->Integral(xmin,xmax);
            }

        }
        std::cout << "x1:    " <<  final_x1_95 << "\tR(x1): " << f2_sys->Eval(final_x1_95) << std::endl;
        std::cout << "x2:    " <<  final_x2_95 << "\tR(x2): " << f2_sys->Eval(final_x2_95) << std::endl;
        std::cout <<"Integral: "<< integral_95<< std::endl;

        std::cout << "x1:    " <<  final_x1_90 << "\tR(x1): " << f2_sys->Eval(final_x1_90) << std::endl;
        std::cout << "x2:    " <<  final_x2_90 << "\tR(x2): " << f2_sys->Eval(final_x2_90) << std::endl;
        std::cout <<"Integral: "<< integral_90<< std::endl;


        TPaveLabel *pav1 = new TPaveLabel(0.5,0.89,0.93,0.94,Form("#bf{Mean: } %f #pm %f", mean,  meanerr),"NDC");
        TPaveLabel *pav2 = new TPaveLabel(0.5,0.81,0.93,0.86,Form("#bf{Sigma:} %f #pm %f",sigma, sigmaerr),"NDC");
        pav1->SetBorderSize(0);   pav2->SetBorderSize(0);
        pav1->SetFillStyle(1001); pav2->SetFillStyle(1001);
        pav1->SetFillColor(0);    pav2->SetFillColor(0); 
        pav1->SetTextFont(12);    pav2->SetTextFont(12);  
        pav1->SetTextSize(0.9);   pav2->SetTextSize(0.9);
        pav1->SetTextAlign(31);   pav2->SetTextAlign(31);
        pav1->SetTextColor(kRed); pav2->SetTextColor(kRed);
        pav1->Draw();             pav2->Draw();

        std::cout << "Mean: " << mean << "\tSigma: " << sigma << "\t 95% limit: "<< mean - 1.64*sigma << std::endl;
        double sigma_sys = 0.788;
        double sigma_total = sqrt(pow(sigma,2) + pow(sigma_sys,2));
        double width_95 = 1.64;
        double width_90 = 1.28;
        //double width_95 = 1.96;
        //double width_90 = 1.64;

        n_gen_line_br[count]     = gen_value;
        //n_fit_line_br[count]     = temp[floor(temp.size()*0.95 + 0.5)-1];
        n_fit_line_br[count]     = mean + width_95*sigma;
        n_fit_line_br[count]     = final_x2_95;
        n_fit_line_br_sys[count] = mean + width_95*sigma_total;
        n_fit_line_br_sys[count] = final_x2_95_sys;

        n_gen_line2_br[count]     = gen_value;
        //n_fit_line2_br[count]     = temp[floor(temp.size()*0.683 + 0.5)-1];
        n_fit_line2_br[count]     = mean + width_90*sigma;
        n_fit_line2_br[count]     = final_x2_90;
        n_fit_line2_br_sys[count] = mean + width_90*sigma_total;
        n_fit_line2_br_sys[count] = final_x2_90_sys;
        
        n_gen_line_br_lower[count]     = gen_value;
        //n_fit_line_br_lower[count]     = temp[floor(temp.size()*0.05 + 0.5)-1];
        n_fit_line_br_lower[count]     = mean - width_95*sigma;
        n_fit_line_br_lower[count]     = final_x1_95;
        n_fit_line_br_sys_lower[count] = mean - width_95*sigma_total;
        n_fit_line_br_sys_lower[count] = final_x1_95_sys;

        n_gen_line2_br_lower[count]     = gen_value;
        //n_fit_line2_br_lower[count]     = temp[floor(temp.size()*0.317 + 0.5)-1]; 
        n_fit_line2_br_lower[count]     = mean - width_90*sigma;  
        n_fit_line2_br_lower[count]     = final_x1_90;  
        n_fit_line2_br_sys_lower[count] = mean - width_90*sigma_total; 
        n_fit_line2_br_sys_lower[count] = final_x1_90_sys; 

        //std::cout << "n_gen_line_br[" << count << "] = " << n_gen_line_br[count] <<std::endl;
        //std::cout << "n_fit_line_br[" << count << "] = " << n_fit_line_br[count] <<std::endl;



        count++;
        std::cout << std::endl;
    
    }


    TCanvas *cav = new TCanvas("Canvas","Sensitivity",200,10,700,500);
  
    TGraph* gr_sensitivity = new TGraph(n_points_new, n_fit_br, n_gen_br);

    gr_sensitivity->SetTitle(Form("Sensitivity %s;Fitted Branching Fraction (#times10^{-7});Generated Branching Fraction (#times10^{-7})",name.c_str())); 
    gr_sensitivity->Draw("AP*");
    
    //TCanvas *cav2 = new TCanvas("Canvas2","Sensitivity2",200,10,700,500);
  
    TGraph* gr_sensitivity_line = new TGraph(n_gen_points, n_fit_line_br, n_gen_line_br);
    gr_sensitivity_line->SetTitle("Sensitivity;Fitted Branching Fraction (#times10^{-7});Generated Branching Fraction (#times10^{-7})"); 
    //gr_sensitivity_line->Draw("AL");  
    
    TGraph* gr_sensitivity_line_sys = new TGraph(n_gen_points, n_fit_line_br_sys, n_gen_line_br);
    gr_sensitivity_line_sys->SetTitle("Sensitivity;Fitted Branching Fraction (#times10^{-7});Generated Branching Fraction (#times10^{-7})"); 
    //gr_sensitivity_line->Draw("AL");  
    
    TGraph* gr_sensitivity_line2 = new TGraph(n_gen_points, n_fit_line2_br, n_gen_line2_br);
    gr_sensitivity_line->SetTitle("Sensitivity;Fitted Branching Fraction (#times10^{-7});Generated Branching Fraction (#times10^{-7})"); 
      
    TGraph* gr_sensitivity_line2_sys = new TGraph(n_gen_points, n_fit_line2_br_sys, n_gen_line2_br);
    gr_sensitivity_line_sys->SetTitle("Sensitivity;Fitted Branching Fraction (#times10^{-7});Generated Branching Fraction (#times10^{-7})"); 
  
    TGraph* gr_sensitivity_line_lower = new TGraph(n_gen_points, n_fit_line_br_lower, n_gen_line_br_lower);
    gr_sensitivity_line_lower->SetTitle("Sensitivity;Fitted Branching Fraction (#times10^{-7});Generated Branching Fraction (#times10^{-7})"); 
    //gr_sensitivity_line->Draw("AL");   
    
    TGraph* gr_sensitivity_line_sys_lower = new TGraph(n_gen_points, n_fit_line_br_sys_lower, n_gen_line_br_lower);
    gr_sensitivity_line_sys_lower->SetTitle("Sensitivity;Fitted Branching Fraction (#times10^{-7});Generated Branching Fraction (#times10^{-7})"); 
    //gr_sensitivity_line->Draw("AL");  
    
    TGraph* gr_sensitivity_line2_lower = new TGraph(n_gen_points, n_fit_line2_br_lower, n_gen_line2_br_lower);
    gr_sensitivity_line2_lower->SetTitle("Sensitivity;Fitted Branching Fraction (#times10^{-7});Generated Branching Fraction (#times10^{-7})"); 
        
    TGraph* gr_sensitivity_line2_sys_lower = new TGraph(n_gen_points, n_fit_line2_br_sys_lower, n_gen_line2_br_lower);
    gr_sensitivity_line2_sys_lower->SetTitle("Sensitivity;Fitted Branching Fraction (#times10^{-7});Generated Branching Fraction (#times10^{-7})"); 
    

    TGraph *grshade = new TGraph(2*n_gen_points);
    for (i=0;i<n_gen_points;i++) {
      grshade->SetPoint(i,n_fit_line_br[i],n_gen_line_br[i]);
      grshade->SetPoint(n_gen_points+i,n_fit_line_br_lower[n_gen_points-i-1],n_gen_line_br_lower[n_gen_points-i-1]);
    }
    grshade->SetFillColor(kYellow);
    TGraph *grshade2 = new TGraph(2*n_gen_points);
    for (i=0;i<n_gen_points;i++) {
      grshade2->SetPoint(i,n_fit_line2_br[i],n_gen_line2_br[i]);
      grshade2->SetPoint(n_gen_points+i,n_fit_line2_br_lower[n_gen_points-i-1],n_gen_line2_br_lower[n_gen_points-i-1]);
    }

    grshade2->SetFillColor(kGreen);

   //grshade->SetFillStyle(3013);
   grshade->SetFillColor(16);
   grshade->Draw("l");
   gr_sensitivity_line->Draw("l");
   gr_sensitivity_line_lower->Draw("l");

    //gr_sensitivity_line->Draw("AL");

    TMultiGraph *mg = new TMultiGraph();
    mg->SetTitle(Form("Sensitivity %s;Fitted #it{B}(%s) (#times10^{-7});Generated #it{B}(%s) (#times10^{-7})",name.c_str(),name.c_str(),name.c_str()));

    mg->Add(grshade,"F");
    mg->Add(grshade2,"F");
    mg->Add(gr_sensitivity_line,"L");
    mg->Add(gr_sensitivity_line_sys,"L");
    mg->Add(gr_sensitivity_line_lower,"L");
    mg->Add(gr_sensitivity_line_sys_lower,"L");

    mg->Add(gr_sensitivity_line2_lower,"L");
    mg->Add(gr_sensitivity_line2,"L");
    //mg->Add(gr_sensitivity_line2_sys_lower,"L");
    //mg->Add(gr_sensitivity_line2_sys,"L");

    
    //mg->Add(gr_sensitivity,"p");
    mg->Draw("A");

    mg->GetXaxis()->SetRangeUser(-4.0,15.0);
    mg->GetYaxis()->SetRangeUser(0.1,7.0);

    gr_sensitivity_line->SetLineColor(kBlack);
    gr_sensitivity_line_sys->SetLineColor(kBlack);
    gr_sensitivity_line_sys->SetLineStyle(kDashed);
    gr_sensitivity_line_lower->SetLineColor(kBlack);
    gr_sensitivity_line_sys_lower->SetLineColor(kBlack);
    gr_sensitivity_line_sys_lower->SetLineStyle(kDashed);

    gr_sensitivity_line2->SetLineColor(kBlack);
    gr_sensitivity_line2_sys->SetLineColor(kBlack);
    gr_sensitivity_line2_sys->SetLineStyle(kDashed);
    gr_sensitivity_line2_lower->SetLineColor(kBlack);
    gr_sensitivity_line2_sys_lower->SetLineColor(kBlack);
    gr_sensitivity_line2_sys_lower->SetLineStyle(kDashed);


    gr_sensitivity->SetMarkerStyle(2);
    //gr_sensitivity_line->SetLineColor(kYellow);
    //gr_sensitivity_line_lower->SetLineColor(kYellow);

    //gr_sensitivity_line2->SetLineColor(kGreen);
    //gr_sensitivity_line2_lower->SetLineColor(kGreen);

    
    
    gr_sensitivity->SetMarkerStyle(2);
    gr_sensitivity_line->SetLineColor(kYellow);
    gr_sensitivity_line->SetFillColor(kYellow);
    gr_sensitivity_line_lower->SetFillColor(0);
    gr_sensitivity_line->Draw("AF");
    gr_sensitivity_line_lower->Draw("F");

    gr_sensitivity_line2->SetLineColor(kGreen);
    gr_sensitivity_line2_lower->SetLineColor(kGreen);



    TLine *line_diagonal = new TLine(0.0,0.0,7.0,7.0); 
    line_diagonal->SetLineStyle(kDashed);
    line_diagonal->SetLineWidth(2);
    line_diagonal->SetLineColor(kWhite);
    line_diagonal->Draw();

    //double limit     = 4.45;
    //double limit_sys = 4.85;

    //double limit     = 3.8;
    //double limit_sys = 4.2;

    double limit     = 4.6;
    double limit_sys = 4.9;

    TLine *line = new TLine(1.16,0.1,1.16,limit_sys);
    line->SetLineWidth(2);
    line->SetLineColor(kRed);
    line->Draw();    

    TLine *line2 = new TLine(-4.0,limit,1.16,limit);
    line2->SetLineStyle(kDashed);
    line2->SetLineWidth(2);
    line2->SetLineColor(kRed);
    line2->Draw();


    TLine *line4 = new TLine(-4.0,limit_sys,1.16,limit_sys);
    line4->SetLineStyle(kDashed);
    line4->SetLineWidth(2);
    line4->SetLineColor(kRed);
    line4->Draw();



    gPad->RedrawAxis();

    TLegend* leg = new TLegend(0.65,0.2,0.92,0.4);
    //TLegend* leg = new TLegend(0.2,0.44,0.5,0.65);
    leg->AddEntry(grshade2,"90% CL Band","f");
    leg->AddEntry(grshade, "95% CL Band","f");
    leg->AddEntry(gr_sensitivity_line_sys,"With Systematics","l");
    leg->Draw();

    cav->Print("results/Sensitivity_plot.eps");
    cav->Print("results/Sensitivity_plot.pdf");

*/
}

std::pair<double,double> GetLimits(double gen_br,
                                   double p0_mean,
                                   double p1_mean,
                                   double p0_stat,
                                   double p1_stat,
                                   double p0_sys,
                                   double p1_sys,
                                   double CL){
    double xmin = -100.0;
    double xmax =  100.0;

    //std::cout << "Calculating limit for CL: " << CL <<std::endl;
    double mean = p0_mean + p1_mean * gen_br; 
    double stat = sqrt(p0_stat + p1_stat * gen_br); 
    double sys  = p0_sys  + p1_sys  * gen_br; 

    double total_sigma = sqrt(pow(stat,2)+pow(sys,2));
    //std::cout << "Using mean:  " << mean <<"\t Stat: "<< stat << "\tSys: " << sys  << "\tTotal: " << total_sigma <<std::endl;
    
    //TCanvas *cav = new TCanvas("Canvas_R","Canvas_R",200,10,700,500);

    //TF1 *f1 = new TF1("f1" , "(x>0)*(pow(x-[0],2)/pow([1],2)) + (x<0)*((pow([0],2) - 2*[0]*x)/pow([1],2))",xmin,xmax);
    TF1 *f2 = new TF1("f2" , "(x>0)*exp(-0.5*(pow(x-[0],2)/pow([1],2))) + (x<0)*exp(-0.5*((pow([0],2) - 2*[0]*x)/pow([1],2)))",xmin,xmax);
    
    f2->SetParameter(0,mean);
    f2->SetParameter(1,total_sigma);
    //f2->Draw();
    double stepsize = 0.0001;
    double maxR = 1.0;
    double minR = 0.0001;
    double integral = -999.0;
    double x1 = -999.0;
    double x2 =  999.0;
    bool found_itergral = false;

    for(double i = 0; i< (int)((maxR-minR)/stepsize); i++){
        double R_temp = i*stepsize+minR;
        double x1_temp = Getx1(R_temp,mean,total_sigma);
        double x2_temp = Getx2(R_temp,mean,total_sigma);
        double integral = f2->Integral(x1_temp,x2_temp)/f2->Integral(xmin,xmax);
        //std::cout << "Checking R = " << R_temp << "\tx1: " << x1_temp<< "\tx2: " << x2_temp << "\tIntegral: "<<  integral <<std::endl;
        if(integral < CL && !found_itergral){
            x1 = x1_temp;
            x2 = x2_temp;
            found_itergral = true;
            //std::cout << "===> Using Values x1: " <<  x1_temp<< "\tx2: " << x2_temp << "\tIntegral: "<<  integral <<std::endl;
            break;

        }
    }

    //delete f1;
    std::pair<double,double> limits;
    limits.first  = x1;
    limits.second = x2;
    return limits;
}

double Getx1(double R, double mean, double sigma){
    double x1_from_gaussian = mean - sqrt(-2*pow(sigma,2)*log(R));

    if(x1_from_gaussian < 0.0){
        return (pow(mean,2) + 2*pow(sigma,2)*log(R));

    } else {
        return x1_from_gaussian;
    }
}

double Getx2(double R, double mean, double sigma){
    return mean + sqrt(-2*pow(sigma,2)*log(R));  
}

void SetLHCbStyle(std::string type = ""){

  if(type!="norm" && type!= "cont" && type!= "surf")
  {
    std::cout <<  "Type must be one of 1D, cont or surf. Defaulting to 1D" <<std::endl;
    type = "norm";
  }

  // Copy and paste from lhcbstyle.C:

  // Use times new roman, precision 2 
  Int_t lhcbFont        = 132;  // Old LHCb style: 62;
  // Line thickness
  Double_t lhcbWidth    = 2.00; // Old LHCb style: 3.00;
  // Text size
  Double_t lhcbTSize    = 0.06; 
  
  // use plain black on white colors
  gROOT->SetStyle("Plain"); 
  TStyle *lhcbStyle= new TStyle("lhcbStyle","LHCb plots style");
  
  //lhcbStyle->SetErrorX(0); //  don't suppress the error bar along X

  lhcbStyle->SetFillColor(1);
  lhcbStyle->SetFillStyle(1001);   // solid
  lhcbStyle->SetFrameFillColor(0);
  lhcbStyle->SetFrameBorderMode(0);
  lhcbStyle->SetPadBorderMode(0);
  lhcbStyle->SetPadColor(0);
  lhcbStyle->SetCanvasBorderMode(0);
  lhcbStyle->SetCanvasColor(0);
  lhcbStyle->SetStatColor(0);
  lhcbStyle->SetLegendBorderSize(0);

  // If you want the usual gradient palette (blue -> red)
  lhcbStyle->SetPalette(1);
  // If you want colors that correspond to gray scale in black and white:
  int colors[8] = {0,5,7,3,6,2,4,1};
  lhcbStyle->SetPalette(8,colors);

  // set the paper & margin sizes
  lhcbStyle->SetPaperSize(20,26);
  lhcbStyle->SetPadTopMargin(0.05);
  lhcbStyle->SetPadRightMargin(0.15); // increase for colz plots
  lhcbStyle->SetPadBottomMargin(0.16);
  lhcbStyle->SetPadLeftMargin(0.14);
  
  // use large fonts
  lhcbStyle->SetTextFont(lhcbFont);
  lhcbStyle->SetTextSize(lhcbTSize);
  lhcbStyle->SetLabelFont(lhcbFont,"x");
  lhcbStyle->SetLabelFont(lhcbFont,"y");
  lhcbStyle->SetLabelFont(lhcbFont,"z");
  lhcbStyle->SetLabelSize(lhcbTSize,"x");
  lhcbStyle->SetLabelSize(lhcbTSize,"y");
  lhcbStyle->SetLabelSize(lhcbTSize,"z");
  lhcbStyle->SetTitleFont(lhcbFont);
  lhcbStyle->SetTitleFont(lhcbFont,"x");
  lhcbStyle->SetTitleFont(lhcbFont,"y");
  lhcbStyle->SetTitleFont(lhcbFont,"z");
  lhcbStyle->SetTitleSize(1.2*lhcbTSize,"x");
  lhcbStyle->SetTitleSize(1.2*lhcbTSize,"y");
  lhcbStyle->SetTitleSize(1.2*lhcbTSize,"z");

  // use medium bold lines and thick markers
  lhcbStyle->SetLineWidth(lhcbWidth);
  lhcbStyle->SetFrameLineWidth(lhcbWidth);
  lhcbStyle->SetHistLineWidth(lhcbWidth);
  lhcbStyle->SetFuncWidth(lhcbWidth);
  lhcbStyle->SetGridWidth(lhcbWidth);
  lhcbStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes
  lhcbStyle->SetMarkerStyle(20);
  lhcbStyle->SetMarkerSize(1.0);

  // label offsets
  lhcbStyle->SetLabelOffset(0.010,"X");
  lhcbStyle->SetLabelOffset(0.010,"Y");

  // by default, do not display histogram decorations:
  lhcbStyle->SetOptStat(0);  
  //lhcbStyle->SetOptStat("emr");  // show only nent -e , mean - m , rms -r
  // full opts at http://root.cern.ch/root/html/TStyle.html#TStyle:SetOptStat
  lhcbStyle->SetStatFormat("6.3g"); // specified as c printf options
  lhcbStyle->SetOptTitle(0);
  lhcbStyle->SetOptFit(0);
  //lhcbStyle->SetOptFit(1011); // order is probability, Chi2, errors, parameters
  //titles
  lhcbStyle->SetTitleOffset(0.95,"X");
  lhcbStyle->SetTitleOffset(0.95,"Y");
  lhcbStyle->SetTitleOffset(1.2,"Z");
  lhcbStyle->SetTitleFillColor(0);
  lhcbStyle->SetTitleStyle(0);
  lhcbStyle->SetTitleBorderSize(0);
  lhcbStyle->SetTitleFont(lhcbFont,"title");
  lhcbStyle->SetTitleX(0.0);
  lhcbStyle->SetTitleY(1.0); 
  lhcbStyle->SetTitleW(1.0);
  lhcbStyle->SetTitleH(0.05);
  
  // look of the statistics box:
  lhcbStyle->SetStatBorderSize(0);
  lhcbStyle->SetStatFont(lhcbFont);
  lhcbStyle->SetStatFontSize(0.05);
  lhcbStyle->SetStatX(0.9);
  lhcbStyle->SetStatY(0.9);
  lhcbStyle->SetStatW(0.25);
  lhcbStyle->SetStatH(0.15);

  // put tick marks on top and RHS of plots
  lhcbStyle->SetPadTickX(1);
  lhcbStyle->SetPadTickY(1);

  // histogram divisions: only 5 in x to avoid label overlaps
  lhcbStyle->SetNdivisions(505,"x");
  lhcbStyle->SetNdivisions(510,"y");
 
   // Extra options
  if(type == "cont"||type=="surf") lhcbStyle->SetPalette(112);
  if(type == "surf") lhcbStyle->SetTitleOffset(1.5,"Y");
  if(type == "surf") lhcbStyle->SetTitleOffset(1.5,"X");
  if(type == "surf") lhcbStyle->SetTitleOffset(0.83,"Z");
  if(type == "norm") lhcbStyle->SetPadRightMargin(0.05);
  

  gROOT->SetStyle("lhcbStyle");
  gROOT->ForceStyle();

  /*
  // add LHCb label
  lhcbName = new TPaveText(gStyle->GetPadLeftMargin() + 0.05,
                           0.87 - gStyle->GetPadTopMargin(),
                           gStyle->GetPadLeftMargin() + 0.20,
                           0.95 - gStyle->GetPadTopMargin(),
                           "BRNDC");
  lhcbName->AddText("LHCb");
  lhcbName->SetFillColor(0);
  lhcbName->SetTextAlign(12);
  lhcbName->SetBorderSize(0);
  */
  TText *lhcbLabel = new TText();
  lhcbLabel->SetTextFont(lhcbFont);
  lhcbLabel->SetTextColor(1);
  lhcbLabel->SetTextSize(lhcbTSize);
  lhcbLabel->SetTextAlign(12);

  TLatex *lhcbLatex = new TLatex();
  lhcbLatex->SetTextFont(lhcbFont);
  lhcbLatex->SetTextColor(1);
  lhcbLatex->SetTextSize(lhcbTSize);
  lhcbLatex->SetTextAlign(12);

  std::cout << "-------------------------" << std::endl;  
  std::cout << "Set LHCb Style - Feb 2012" << std::endl;  
  std::cout << "Mode : " << type           << std::endl;
  std::cout << "-------------------------" << std::endl;  


}
