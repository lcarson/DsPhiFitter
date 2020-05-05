void SetLHCbStyle(std::string);
#include "tools.h"
#include <stdlib.h> 
// ===========================================================================

void plotSensitivity(std::string mode = "DsPhi"){


    std::string name;
    if (mode == "DsPhi") name = "#it{B^{+}#rightarrowD_{s}^{+}#phi}";
    if (mode == "DKst0") name = "it{B^{+}#rightarrowD^{+}K^{*0}}";

    std::cout << "Using mode name: " << mode << " Name: " << name << std::endl;

    std::string sensitivityDir = "sensitivityDir/";
    std::vector<std::string> fileName;
    getdir(sensitivityDir,fileName);
    std::map<double,int> gen_count;
    std::map<double,std::vector<double>> gen_br;

    int n_points = 0;
    for(std::vector<std::string>::iterator it=fileName.begin();it!=fileName.end();it++){
        if( (*it).find("sensitivity")  ==std::string::npos ) continue;
        n_points++;
    }

    double* n_gen_br   = new double[n_points];
    double* n_fit_br   = new double[n_points];
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
        double fitted_br = param_final->getVal();
        if(count%1000==0)  std::cout << "Br: " << br << " Fitted Br: " << fitted_br << " seed: " << seed << " Number: " << num << std::endl;
        n_gen_br[count] = br;
        n_fit_br[count] = fitted_br;
        count++ ;
        n_points_new++;
        gen_count[br]++;
        gen_br[br].push_back(fitted_br);


        /*

        TIter nextkey(newfile.GetListOfKeys());
        TKey *key;
        while((key=(TKey*)nextkey())){
          if(std::string(key->GetClassName())!=std::string("RooFitResult")) continue;
          RooFitResult *r = (RooFitResult*)key->ReadObj();
          //if(std::string(r->GetName()) != std::string("toy_withSig")) continue;
          //if( !newfile->GetListOfKeys()->Contains("toy_withSig"))  ) continue; 
          std::cout << "RooFitResult Name: " << std::endl;
          r->Print();
          nFits++;      
          if(r->covQual()<2){nBad++;continue;}
          if(r->covQual()<3){nFPD++;continue;}
          if(r->status()!=0){nBadMINOS++;}//continue;}
          
          RooRealVar* param_final = (RooRealVar*) (r->floatParsFinal()).find("Branching_fraction");
          double fitted_br = param_final->getVal();

        

        
        TFile *newfile = TFile::Open((sensitivityDir+*it).c_str());
        if(!newfile){
          std::cout << "File not found: " << sensitivityDir+*it << std::endl;
          continue;
        }
        RooFitResult* result = (RooFitResult*)newfile->Get("toy_withSig");
        RooRealVar* param_final = (RooRealVar*) (result->floatParsFinal()).find("Branching_fraction");
        double fitted_br = param_final->getVal();
        

          std::cout << "Br: " << br << " Fitted Br: " << fitted_br << " seed: " << seed << " Number: " << num << std::endl;
          n_gen_br[count] = br;
          n_fit_br[count] = fitted_br;
          count++ ;
          gen_count[br]++;
          gen_br[br].push_back(fitted_br);
          //newfile->Close();
        }
      */

    }

    int n_gen_points = gen_br.size();
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


    count = 0;
    SetLHCbStyle("norm");


    for(std::map<double,std::vector<double>>::iterator it=gen_br.begin();it!=gen_br.end();it++){ 
        std::vector<double> temp = it->second;
        double gen_value = it->first;
        std::sort (temp.begin(), temp.end());

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


        TPaveLabel *pav1 = new TPaveLabel(0.2,0.89,0.93,0.94,Form("#bf{Mean: } %f #pm %f", mean,  meanerr),"NDC");
        TPaveLabel *pav2 = new TPaveLabel(0.2,0.81,0.93,0.86,Form("#bf{Sigma:} %f #pm %f",sigma, sigmaerr),"NDC");
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
/*
   //grshade->SetFillStyle(3013);
   grshade->SetFillColor(16);
   grshade->Draw("l");
   gr_sensitivity_line->Draw("l");
   gr_sensitivity_line_lower->Draw("l");
*/
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
    mg->GetYaxis()->SetRangeUser(0.2,7.0);
    mg->GetXaxis()->SetTitleSize(0.05);
    mg->GetYaxis()->SetTitleSize(0.05);
    mg->GetXaxis()->SetTitleOffset(1.2);
    mg->GetYaxis()->SetTitleOffset(1.2);
    mg->GetXaxis()->SetLabelSize(0.05);
    mg->GetYaxis()->SetLabelSize(0.05);

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

    /*
    
    gr_sensitivity->SetMarkerStyle(2);
    gr_sensitivity_line->SetLineColor(kYellow);
    gr_sensitivity_line->SetFillColor(kYellow);
    gr_sensitivity_line_lower->SetFillColor(0);
    gr_sensitivity_line->Draw("AF");
    gr_sensitivity_line_lower->Draw("F");

    gr_sensitivity_line2->SetLineColor(kGreen);
    gr_sensitivity_line2_lower->SetLineColor(kGreen);

*/

    TLine *line_diagonal = new TLine(0.0,0.0,7.0,7.0); 
    line_diagonal->SetLineStyle(kDashed);
    line_diagonal->SetLineWidth(2);
    line_diagonal->SetLineColor(kWhite);
    line_diagonal->Draw();

    //double limit     = 4.45;
    //double limit_sys = 4.85;

    //double limit     = 3.8;
    //double limit_sys = 4.2;

    double limit     = 4.4;
    double limit_sys = 4.9;

    TLine *line = new TLine(1.16,0.2,1.16,limit_sys);
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

    TLegend* leg = new TLegend(0.58,0.2,0.92,0.4);
    //TLegend* leg = new TLegend(0.2,0.44,0.5,0.65);
    leg->AddEntry(grshade2,"90% CL band","f");
    leg->AddEntry(grshade, "95% CL band","f");
    leg->AddEntry(gr_sensitivity_line_sys,"With syst. uncertainty","l");
    leg->Draw();

    cav->Print("results/Sensitivity_plot.eps");
    cav->Print("results/Sensitivity_plot.pdf");
    cav->Print("results/Sensitivity_plot.C");
    cav->Print("results/Sensitivity_plot.root");
    cav->Print("results/Sensitivity_plot.png");
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
