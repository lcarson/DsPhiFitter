void SetLHCbStyle(std::string);
#include "tools2.h"
#include <stdlib.h> 
// ===========================================================================

void plotToys(){
    SetLHCbStyle("oneD");
    


    std::map<std::string,std::string> proper_names;
    proper_names["yield_peak_DsD0_Ds2PhiPi_toy_both_DsBDTbin1_PhiBDTbin1_both_both"]            = "N(#it{B^{+}#rightarrow}(#it{D_{s}^{+}#rightarrow#phi#pi^{+}})#it{#bar{D}^{0}})";       
    proper_names["yield_peak_DsD0_Ds2KKPi_toy_both_DsBDTbin1_PhiBDTbin1_both_both"]             = "N(#it{B^{+}#rightarrow}(#it{D_{s}^{+}#rightarrowK^{+}K^{#minus}#pi^{+}})#it{#bar{D}^{0}})";        
    proper_names["yield_peak_DsD0_Ds2PiPiPi_toy_both_DsBDTbin1_PhiBDTbin1_both_both"]           = "N(#it{B^{+}#rightarrow}(#it{D_{s}^{+}#rightarrow#pi^{+}#pi^{#minus}#pi^{+}})#it{#bar{D}^{0}})";       
    proper_names["yield_peak_DsD0_Ds2KPiPi_toy_both_DsBDTbin1_PhiBDTbin1_both_both"]            = "N(#it{B^{+}#rightarrow}(#it{D_{s}^{+}#rightarrowK^{+}#pi^{#minus}#pi^{+}})#it{#bar{D}^{0}})";       
    proper_names["yield_peak_total_DsPhi_Ds2PhiPi_toy_both_DsBDTbin1_PhiBDTbin1_both_both"]     = "N(#it{B^{+}#rightarrow}(#it{D_{s}^{+}#rightarrow#phi#pi^{+}})#it{#phi})";             
    proper_names["yield_peak_total_DsPhi_Ds2KKPi_toy_both_DsBDTbin1_PhiBDTbin1_both_both"]      = "N(#it{B^{+}#rightarrow}(#it{D_{s}^{+}#rightarrowK^{+}K^{#minus}#pi^{+}})#it{#phi})";              
    proper_names["yield_peak_total_DsPhi_Ds2PiPiPi_toy_both_DsBDTbin1_PhiBDTbin1_both_both"]    = "N(#it{B^{+}#rightarrow}(#it{D_{s}^{+}#rightarrow#pi^{+}#pi^{#minus}#pi^{+}})#it{#phi})";               
    proper_names["yield_peak_total_DsPhi_Ds2KPiPi_toy_both_DsBDTbin1_PhiBDTbin1_both_both"]     = "N(#it{B^{+}#rightarrow}(#it{D_{s}^{+}#rightarrowK^{+}#pi^{#minus}#pi^{+}})#it{#phi})";                 
    proper_names["Branching_fraction"] = "B(#it{B^{+}#rightarrowD_{s}^{+}#phi}) (x 10^{-7})";
    std::string dir = "toysDir/";
    
    std::vector<std::string> files;
    getdir(dir,files);
    std::cout << "Number of files: " << files.size() << std::endl; 
    if (files.size()==0) return;
    
    int nFits=0;
    int nBad=0;
    int nBadMINOS=0;
    int nGood=0;
    int nFPD=0;


    std::vector<RooFitResult> Results_DsKK; 

    std::cout << "Looking for DsPhi Toys..." <<std::endl;
    int count =0;
    for(std::vector<std::string>::iterator it=files.begin();it!=files.end();it++){ 
        if( (*it).find("toy")  ==std::string::npos ) continue;
        //std::cout << "Processing file : " << (*it) << std::endl;
        std::string filename = (*it);
        
        TFile f((dir+filename).c_str(),"READ");
        TIter nextkey(f.GetListOfKeys());
        TKey *key;
        while((key=(TKey*)nextkey())){
            if(std::string(key->GetClassName())!=std::string("RooFitResult")) continue;
            RooFitResult *r = (RooFitResult*)key->ReadObj();
            nFits++;      
            if(r->covQual()<2){nBad++;continue;}
            if(r->covQual()<3){nFPD++;continue;}
            if(r->status()!=0){nBadMINOS++;}
            nGood++;
            Results_DsKK.push_back(*r);
        }    
        std::cout <<"\r   Collecting "<<nFits<<" toys, of which "<<nBad<<" don't converge and "<<100*float(nFPD)/nFits<<"% are forced positive definite and "<<100*float(nBadMINOS)/nFits<<"% have MINOS problems. "<< nGood << " are good."<<std::flush;

    }   
    std::cout<<std::endl;
    if (Results_DsKK.size()==0) return;
    int total = Results_DsKK.size();


    // Get floating vars 

    std::map<std::string,TH1D*> h_Value_DsKK;
    std::map<std::string,TH1D*> h_Error_DsKK;
    std::map<std::string,TH1D*> h_Error_hi_DsKK;
    std::map<std::string,TH1D*> h_Error_lo_DsKK;
    std::map<std::string,TH1D*> h_Pull_DsKK;
    std::map<std::string,double> Init_DsKK;

    // Make histograms
    RooArgList finalPars =  Results_DsKK[0].floatParsFinal();

    RooRealVar* temp_var=0;
    TIterator* it = finalPars.createIterator();
    while((temp_var = (RooRealVar*)it->Next())){
         std::string name = temp_var->GetName();
        
        std::cout << "Final: "      << temp_var->getVal();
        std::cout << "\t+/- "       << temp_var->getError();
        std::cout << "\t Min: "     << temp_var->getMin();
        std::cout << "\t Max: "     << temp_var->getMax();
        std::cout << "\tVariable: " << temp_var->GetName();
        //std::cout <<std::endl;

        h_Value_DsKK[name] = new TH1D(Form("h_%s_Value",name.c_str()), Form("%s Value",name.c_str()),20,   0,  0    );
        h_Error_DsKK[name] = new TH1D(Form("h_%s_Error",name.c_str()), Form("%s Error",name.c_str()),20,   0,  0    );        
        h_Pull_DsKK[name]  = new TH1D(Form("h_%s_Pull",name.c_str()),  Form("%s Pull",name.c_str()), 20,-5.0,5.0    );

        h_Error_hi_DsKK[name] = new TH1D(Form("h_%s_Error_hi",name.c_str()), Form("%s Error_hi",name.c_str()),20,   0,  0    );        
        h_Error_lo_DsKK[name] = new TH1D(Form("h_%s_Error_lo",name.c_str()), Form("%s Error_lo",name.c_str()),20,   0,  0    );        

        std::string text_filename = Form("%s/gen_vals/Init_Val_%s.txt",dir.c_str(),name.c_str());
        //std::cout  << "File name:" << text_filename << std::endl;
        std::ifstream input;
        input.open(text_filename.c_str(),std::ifstream::in);
        if(!input){
            std::cout<<" Can't find file " << text_filename << std::endl; return;
        }
        std::string line = "";
        getline(input,line);
        Init_DsKK[name] = atof(line.c_str());
        std::cout  << "\t Init: " << Init_DsKK[name] << std::endl;
        input.close();

    }

    // Calculate pulls and fill histogram
    int i = 0;
    for(std::vector<RooFitResult>::iterator it_result=Results_DsKK.begin();it_result!=Results_DsKK.end();it_result++){
        if(i%100==0) std::cout << "Progress: " << 100*(double)i/(double)total << "%"<<  std::endl;  
        i++;
        RooArgList finalPars_loop =  (*it_result).floatParsFinal();

        RooRealVar* temp_var=0;
        TIterator* it_loop = finalPars_loop.createIterator();
        while((temp_var = (RooRealVar*)it_loop->Next())){
            std::string name = temp_var->GetName();
            double value    = temp_var->getVal();
            double error    = temp_var->getError();
            double error_hi = temp_var->getAsymErrorHi();
            double error_lo = temp_var->getAsymErrorLo();


            double pull = (temp_var->getVal()-Init_DsKK[name])/temp_var->getError();
            double asym_pull = 0;
            if(value < Init_DsKK[name]) {
                asym_pull = (value -Init_DsKK[name])/error_hi ;
            }else {
                asym_pull = -(value - Init_DsKK[name])/error_lo ;
            }

            if(name == "Branching_fraction" && error > 1000.0) continue;
            h_Value_DsKK[name]->Fill(value);
            h_Error_DsKK[name]->Fill(error);
            h_Error_hi_DsKK[name]->Fill(error_hi);
            h_Error_lo_DsKK[name]->Fill(error_lo);
            h_Pull_DsKK[name]->Fill(asym_pull);
        }
    }


    // Plot histograms
    temp_var=0;
    it = finalPars.createIterator();
    while((temp_var = (RooRealVar*)it->Next())){
        std::string name = temp_var->GetName();

        //TCanvas* can = new TCanvas(Form("can_%s",name.c_str()),Form("can_%s",name.c_str()), 900, 300);
        TCanvas* can = new TCanvas(Form("can_%s_Value",name.c_str()),Form("can_%s",name.c_str()));
        //can->Divide(3,1);
        
        //can->cd(1);
        h_Value_DsKK[name]->Draw("PE0");
        //h_Value_DsKK[name]->SetTitle(";Value;Number of Toys");
        h_Value_DsKK[name]->SetTitle(Form(";%s Value;Number of Toys",proper_names[name].c_str()));

        //can->cd(2);
        TCanvas* can2 = new TCanvas(Form("can_%s_Error",name.c_str()),Form("can_%s",name.c_str()));
        h_Error_DsKK[name]->Draw("PE0"); 
        h_Error_DsKK[name]->SetTitle(Form(";%s Error;Number of Toys",proper_names[name].c_str())); 
        
        can->cd(3);
        TCanvas* can3 = new TCanvas(Form("can_%s_Pull",name.c_str()),Form("can_%s",name.c_str()));
        h_Pull_DsKK[name]->Draw("PE0");
        h_Pull_DsKK[name]->SetTitle(Form(";%s Pull;Number of Toys",proper_names[name].c_str())); 


        TF1 *gaussian = new TF1("Gaussian","gaus",-5.0,5.0);
        h_Pull_DsKK[name]->Fit("Gaussian","Q0R");
        double mean  = gaussian->GetParameter(1);
        double sigma = gaussian->GetParameter(2);
        const double* err=gaussian->GetParErrors();
        double meanerr=err[1];
        double sigmaerr=err[2];
        //can->cd(3) ; 
        gaussian->SetLineColor(kBlue); 
        gaussian->Draw("same");

        TPaveLabel *pav1 = new TPaveLabel(0.2,0.89,0.93,0.94,Form("#bf{Mean: } %.3f #pm %.3f", mean,  meanerr),"NDC");
        TPaveLabel *pav2 = new TPaveLabel(0.2,0.81,0.93,0.86,Form("#bf{Sigma:} %.3f #pm %.3f",sigma, sigmaerr),"NDC");
        
        pav1->SetBorderSize(0);   pav2->SetBorderSize(0);
        pav1->SetFillStyle(1001); pav2->SetFillStyle(1001);
        pav1->SetFillColor(0);    pav2->SetFillColor(0); 
        pav1->SetTextFont(12);    pav2->SetTextFont(12);  
        pav1->SetTextSize(0.9);   pav2->SetTextSize(0.9);
        pav1->SetTextAlign(31);   pav2->SetTextAlign(31);
        pav1->SetTextColor(kRed); pav2->SetTextColor(kRed);
        
        //pav1->Draw();             pav2->Draw();

        TLatex *t1 = new TLatex(0.92,0.85,Form("#bf{Mean:} %.3f #pm %.3f", mean,  meanerr));
        TLatex *t2 = new TLatex(0.92,0.75,Form("#bf{Sigma:} %.3f #pm %.3f",sigma, sigmaerr));
        t1->SetNDC();             t2->SetNDC();
        t1->SetTextColor(kRed);   t2->SetTextColor(kRed);
        //t1->SetTextSize(0.9);   t2->SetTextSize(0.9);
        t1->SetTextAlign(31);   t2->SetTextAlign(31);
        t1->Draw();             t2->Draw();



        can->Print(Form("toysDir/plots/Plots_DsKK_Value_%s.eps",name.c_str()));
        can->Print(Form("toysDir/plots/Plots_DsKK_Value_%s.pdf",name.c_str()));

        can2->Print(Form("toysDir/plots/Plots_DsKK_Error_%s.eps",name.c_str()));
        can2->Print(Form("toysDir/plots/Plots_DsKK_Error_%s.pdf",name.c_str()));

        can3->Print(Form("toysDir/plots/Plots_DsKK_Pull_%s.eps",name.c_str()));
        can3->Print(Form("toysDir/plots/Plots_DsKK_Pull_%s.pdf",name.c_str()));

    }

    std::cout << "Done" <<std::endl; 
}

