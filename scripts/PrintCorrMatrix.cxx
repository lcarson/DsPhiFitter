void PrintCorrMatrix(std::string filelocation = "/home/hadavizadeh/DST/B2DsPhi/Systematics/2017-03-14_Nominal_fit_Merged_DsPhi/Bu2DsPhi_Fitter_copy/systematicsDir/vary_lastRooFitResult_summed_splitHel_s21_s21r1_s24_s26_Ds2PhiPi_Ds2KKPi_Ds2PiPiPi_Ds2KPiPi_SEED_700_16243.root"){
    
    //filelocation = "/home/hadavizadeh/DST/B2DsPhi/Systematics/2017-04-25_fix_DsKK_to_DsD0_DsPhi/Bu2DsPhi_Fitter_copy/systematicsDir/vary_lastRooFitResult_summed_splitHel_s21_s21r1_s24_s26_Ds2PhiPi_Ds2KKPi_Ds2PiPiPi_Ds2KPiPi_SEED_4000_13705.root";
    filelocation = "/home/hadavizadeh/data/B2DsPhi/Systematics/2017-05-12_fully_unblind_DsPhi/Bu2DsPhi_Fitter_copy/systematicsDir/vary_lastRooFitResult_summed_splitHel_s21_s21r1_s24_s26_Ds2PhiPi_Ds2KKPi_Ds2PiPiPi_Ds2KPiPi_SEED_5000_13704.root";
    
    bool limited = true;

    double threshold = 0.2;

    std::map<std::string,std::string> Description;
    // ============== >>
    Description["Branching_fraction"]           = "Branching fraction (#times10^{-7})";
    Description["DsKK_to_DsD0_Ratio"]           = "Ratio of D_{s}K^{+}K^{-} events to D_{s} D^{0}";
    Description["Dsa1_to_DsstPhi_fraction"]     = "Fraction of D_{s}^{*} #phi in (D_{s}^{(*)}KK^{*0} + D_{s}^{*} #phi)";
    Description["Fraction_DsstD0_Helbin1"]      = "Fraction of D_{s}^{*}D^{0} in (D_{s}^{*}D^{0}+D_{s}D^{*0}) in 1st helicity bin  ";
    Description["Fraction_DsstD0_Helbin2"]      = "Fraction of D_{s}^{*}D^{0} in (D_{s}^{*}D^{0}+D_{s}D^{*0}) in 2nd helicity bin  ";
    Description["Ratio_DD_to_Dsa1"]             = "Ratio of DD^{'} backgrounds to (D_{s}^{*} #phi + D_{s}^{(*)}KK^{*0})";
    Description["Ratio_DsstDst0_to_Lowmass"]    = "Ratio of D_{s}^{*}D^{*0} to (D_{s}^{*}D^{0} + D_{s}D^{*0})";
    Description["global_comb_slope"]            = "Combinatorial slope";
    Description["frac_Dsa1_DsstKKst"]           = "Fraction of D_{s}KK^{*0} in (D_{s}^{*}KK^{*0}+D_{s}KK^{*0})";
    Description["global_csi"]                   = "Relative heights double peaked part reco shapes";
    Description["global_mean"]                  = "Mean B mass (MeV)";
    Description["global_shift"]                 = "Part reco background mass offset (MeV)";
    Description["low_mass_total_DsD0_Ds2KKPi_both_both_DsBDTbin1_PhiBDTbin1_both_both"]    = "Yield of (D_{s}^{*}D^{0} + D_{s}D^{*0}) in D_{s} #rightarrow KK#pi";
    Description["low_mass_total_DsD0_Ds2KPiPi_both_both_DsBDTbin1_PhiBDTbin1_both_both"]   = "Yield of (D_{s}^{*}D^{0} + D_{s}D^{*0}) in D_{s} #rightarrow K#pi#pi";
    Description["low_mass_total_DsD0_Ds2PhiPi_both_both_DsBDTbin1_PhiBDTbin1_both_both"]   = "Yield of (D_{s}^{*}D^{0} + D_{s}D^{*0}) in D_{s} #rightarrow #phi#pi";
    Description["low_mass_total_DsD0_Ds2PiPiPi_both_both_DsBDTbin1_PhiBDTbin1_both_both"]  = "Yield of (D_{s}^{*}D^{0} + D_{s}D^{*0}) in D_{s} #rightarrow #pi#pi#pi";
    Description["low_mass_total_DsPhi_Ds2KKPi_both_both_DsBDTbin1_PhiBDTbin1_both_both"]   = "Yield of (D_{s}^{*} #phi + D_{s}^{(*)}KK^{*0}) in D_{s} #rightarrow KK#pi";
    Description["low_mass_total_DsPhi_Ds2KPiPi_both_both_DsBDTbin1_PhiBDTbin1_both_both"]  = "Yield of (D_{s}^{*} #phi + D_{s}^{(*)}KK^{*0}) in D_{s} #rightarrow K#pi#pi";
    Description["low_mass_total_DsPhi_Ds2PhiPi_both_both_DsBDTbin1_PhiBDTbin1_both_both"]  = "Yield of (D_{s}^{*} #phi + D_{s}^{(*)}KK^{*0}) in D_{s} #rightarrow #phi#pi";
    Description["low_mass_total_DsPhi_Ds2PiPiPi_both_both_DsBDTbin1_PhiBDTbin1_both_both"] = "Yield of (D_{s}^{*} #phi + D_{s}^{(*)}KK^{*0}) in D_{s} #rightarrow K#pi#pi";
    Description["sigma_DsD0_Ds2KKPi"]           = "#sigma_{1} for D_{s} #rightarrow KK#pi  (MeV)";
    Description["sigma_DsD0_Ds2KPiPi"]          = "#sigma_{1} for D_{s} #rightarrow K#pi#pi (MeV)";
    Description["sigma_DsD0_Ds2PhiPi"]          = "#sigma_{1} for D_{s} #rightarrow #phi#pi  (MeV)";
    Description["sigma_DsD0_Ds2PiPiPi"]         = "#sigma_{1} for D_{s} #rightarrow #pi#pi#pi  (MeV)";
    Description["splitHel_DsD0_PR_peak_fraction"] = "Fraction of normalisation part reco in 1st helicity bin";
    Description["splitHel_DsD0_peak_fraction"]    = "Fraction of normalisation peak in 1st helicity bin";
    
    Description["yield_comb_DsD0_Ds2KKPi_both_Helbin1_DsBDTbin1_PhiBDTbin1_both_both"]        = "Comb. Yield D_{s}D^{0} D_{s} #rightarrow KK#pi in 1st helicity bin";
    Description["yield_comb_DsD0_Ds2KKPi_both_Helbin2_DsBDTbin1_PhiBDTbin1_both_both"]        = "Comb. Yield D_{s}D^{0} D_{s} #rightarrow KK#pi in 2nd helicity bin";
    Description["yield_comb_DsD0_Ds2KPiPi_both_Helbin1_DsBDTbin1_PhiBDTbin1_both_both"]       = "Comb. Yield D_{s}D^{0} D_{s} #rightarrow K#pi#pi in 1st helicity bin";
    Description["yield_comb_DsD0_Ds2KPiPi_both_Helbin2_DsBDTbin1_PhiBDTbin1_both_both"]       = "Comb. Yield D_{s}D^{0} D_{s} #rightarrow K#pi#pi in 2nd helicity bin";
    Description["yield_comb_DsD0_Ds2PhiPi_both_Helbin1_DsBDTbin1_PhiBDTbin1_both_both"]       = "Comb. Yield D_{s}D^{0} D_{s} #rightarrow #phi#pi in 1st helicity bin";
    Description["yield_comb_DsD0_Ds2PhiPi_both_Helbin2_DsBDTbin1_PhiBDTbin1_both_both"]       = "Comb. Yield D_{s}D^{0} D_{s} #rightarrow #phi#pi in 2nd helicity bin";
    Description["yield_comb_DsD0_Ds2PiPiPi_both_Helbin1_DsBDTbin1_PhiBDTbin1_both_both"]      = "Comb. Yield D_{s}D^{0} D_{s} #rightarrow #pi#pi#pi in 1st helicity bin";
    Description["yield_comb_DsD0_Ds2PiPiPi_both_Helbin2_DsBDTbin1_PhiBDTbin1_both_both"]      = "Comb. Yield D_{s}D^{0} D_{s} #rightarrow #pi#pi#pi in 2nd helicity bin";
    
    Description["yield_comb_DsPhiSide_Ds2KKPi_both_Helbin1_DsBDTbin1_PhiBDTbin1_both_both"]   = "Comb. Yield D_{s}#phi sideband D_{s} #rightarrow KK#pi in 1st helicity bin";
    Description["yield_comb_DsPhiSide_Ds2KKPi_both_Helbin2_DsBDTbin1_PhiBDTbin1_both_both"]   = "Comb. Yield D_{s}#phi sideband D_{s} #rightarrow KK#pi in 2nd helicity bin";
    Description["yield_comb_DsPhiSide_Ds2KPiPi_both_Helbin1_DsBDTbin1_PhiBDTbin1_both_both"]  = "Comb. Yield D_{s}#phi sideband D_{s} #rightarrow K#pi#pi in 1st helicity bin";
    Description["yield_comb_DsPhiSide_Ds2KPiPi_both_Helbin2_DsBDTbin1_PhiBDTbin1_both_both"]  = "Comb. Yield D_{s}#phi sideband D_{s} #rightarrow K#pi#pi in 2nd helicity bin";
    Description["yield_comb_DsPhiSide_Ds2PhiPi_both_Helbin1_DsBDTbin1_PhiBDTbin1_both_both"]  = "Comb. Yield D_{s}#phi sideband D_{s} #rightarrow #phi#pi in 1st helicity bin";
    Description["yield_comb_DsPhiSide_Ds2PhiPi_both_Helbin2_DsBDTbin1_PhiBDTbin1_both_both"]  = "Comb. Yield D_{s}#phi sideband D_{s} #rightarrow #phi#pi in 2nd helicity bin";
    Description["yield_comb_DsPhiSide_Ds2PiPiPi_both_Helbin1_DsBDTbin1_PhiBDTbin1_both_both"] = "Comb. Yield D_{s}#phi sideband D_{s} #rightarrow #pi#pi#pi in 1st helicity bin";
    Description["yield_comb_DsPhiSide_Ds2PiPiPi_both_Helbin2_DsBDTbin1_PhiBDTbin1_both_both"] = "Comb. Yield D_{s}#phi sideband D_{s} #rightarrow #pi#pi#pi in 2nd helicity bin";
    
    Description["yield_comb_DsPhi_Ds2KKPi_both_Helbin1_DsBDTbin1_PhiBDTbin1_both_both"]       = "Comb. Yield D_{s}#phi D_{s} #rightarrow KK#pi in 1st helicity bin";
    Description["yield_comb_DsPhi_Ds2KKPi_both_Helbin2_DsBDTbin1_PhiBDTbin1_both_both"]       = "Comb. Yield D_{s}#phi D_{s} #rightarrow KK#pi in 2nd helicity bin";
    Description["yield_comb_DsPhi_Ds2KPiPi_both_Helbin1_DsBDTbin1_PhiBDTbin1_both_both"]      = "Comb. Yield D_{s}#phi D_{s} #rightarrow K#pi#pi in 1st helicity bin";
    Description["yield_comb_DsPhi_Ds2KPiPi_both_Helbin2_DsBDTbin1_PhiBDTbin1_both_both"]      = "Comb. Yield D_{s}#phi D_{s} #rightarrow K#pi#pi in 2nd helicity bin";
    Description["yield_comb_DsPhi_Ds2PhiPi_both_Helbin1_DsBDTbin1_PhiBDTbin1_both_both"]      = "Comb. Yield D_{s}#phi D_{s} #rightarrow #phi#pi in 1st helicity bin";
    Description["yield_comb_DsPhi_Ds2PhiPi_both_Helbin2_DsBDTbin1_PhiBDTbin1_both_both"]      = "Comb. Yield D_{s}#phi D_{s} #rightarrow #phi#pi in 2nd helicity bin";
    Description["yield_comb_DsPhi_Ds2PiPiPi_both_Helbin1_DsBDTbin1_PhiBDTbin1_both_both"]     = "Comb. Yield D_{s}#phi D_{s} #rightarrow #pi#pi#pi in 1st helicity bin";
    Description["yield_comb_DsPhi_Ds2PiPiPi_both_Helbin2_DsBDTbin1_PhiBDTbin1_both_both"]     = "Comb. Yield D_{s}#phi D_{s} #rightarrow #pi#pi#pi in 2nd helicity bin";
    
    Description["yield_peak_DsD0_Ds2KKPi_both_both_DsBDTbin1_PhiBDTbin1_both_both"]           = "Total Normalisation Yield D_{s}D^{0} D_{s} #rightarrow KK#pi";
    Description["yield_peak_DsD0_Ds2KPiPi_both_both_DsBDTbin1_PhiBDTbin1_both_both"]          = "Total Normalisation Yield D_{s}D^{0} D_{s} #rightarrow K#pi#pi";
    Description["yield_peak_DsD0_Ds2PhiPi_both_both_DsBDTbin1_PhiBDTbin1_both_both"]          = "Total Normalisation Yield D_{s}D^{0} D_{s} #rightarrow #phi#pi";
    Description["yield_peak_DsD0_Ds2PiPiPi_both_both_DsBDTbin1_PhiBDTbin1_both_both"]         = "Total Normalisation Yield D_{s}D^{0} D_{s} #rightarrow #pi#pi#pi";



    TFile* inputfile = TFile::Open(filelocation.c_str(),"READ");
    if(!inputfile) {
        std::cout << "Can't find file "<< filelocation << std::endl;
        return;
    }
    std::string fitname = "fitresult_model_reducedData";
    RooFitResult *result = (RooFitResult*) inputfile->Get(fitname.c_str());
    if(!result) {
        std::cout << "Can't find RooFitResult "<< fitname << std::endl;
        return;
    }   

    gStyle->SetOptStat(0);
    TCanvas* corrCanv = new TCanvas("corrCanv","",1000,600); 
    corrCanv->cd();

    if(!limited){
        TH2* corr_martrix = result->correlationHist();
        corr_martrix->Draw("colz");
        corr_martrix->SetMarkerSize(0.8);
        corr_martrix->Draw("text same");
    } else {

        RooArgList arglist = result->floatParsFinal();

        RooRealVar* temp_var1=0;
        TIterator* it = arglist.createIterator();
        std::vector<std::string> selected_vars;
        while((temp_var1 = (RooRealVar*)it->Next())) {
            if(temp_var1->isConstant() or temp_var1->InheritsFrom("RooAbsCategory") ) continue;
            std::string varname1 = temp_var1->GetName();
            std::cout << "Looking at variable: " << varname1 << std::endl;
            bool Correlated = false;
            RooRealVar* temp_var2=0;
            TIterator* it2 = arglist.createIterator();
            while((temp_var2 = (RooRealVar*)it2->Next())) {
              if(temp_var2->isConstant() or temp_var2->InheritsFrom("RooAbsCategory") ) continue;
              
              std::string varname2 = temp_var2->GetName();
              double corr = result->correlation(varname1.c_str(),varname2.c_str());
              //std::cout << "      Comparing to: " << varname2 << " Correlation: " << corr <<std::endl;
              if( (varname1!=varname2) && abs(corr) > threshold) Correlated = true;

            }

            if(Correlated){ 
                std::cout << "  At least one correlation greater than " << threshold << std::endl;
                selected_vars.push_back(varname1);
            } else{
                std::cout << "  No correlations greater than " << threshold << std::endl;
            }

        }
        
        std::cout << "Drawing histogram " << std::endl;

        int n = selected_vars.size();
        TH2D* hh = new TH2D("h2_name","h2_name",n,0,n,n,0,n) ;
        
        for (Int_t i = 0 ; i<n ; i++) {
            for (Int_t j = 0 ; j<n; j++) {
                double temp_corr = result->correlation(selected_vars[i].c_str(),selected_vars[j].c_str());
                hh->Fill(i+0.5,n-j-0.5,temp_corr) ;
            }
            hh->GetXaxis()->SetBinLabel(i+1,Description[selected_vars[i]].c_str()) ;
            hh->GetYaxis()->SetBinLabel(n-i,Description[selected_vars[i]].c_str()) ; 
        }
        
        std::cout << "Done making histogram " << std::endl;

        gStyle->SetPadLeftMargin(0.5);
        gStyle->SetPadBottomMargin(0.5);
        hh->SetMinimum(-1) ;
        hh->SetMaximum(+1) ;
        gPad->SetTopMargin( 0.05);
        gPad->SetLeftMargin( 0.15);
        gPad->SetBottomMargin(0.1);
        gPad->SetRightMargin(0.15);
        hh->GetXaxis()->RotateTitle(false);
        hh->GetXaxis()->SetLabelSize(0.01);
        hh->GetYaxis()->SetLabelSize(0.015);
        hh->SetTitle("");
        hh->SetMarkerSize(0.5);
        gStyle->SetPaintTextFormat("4.1f");
        //gStyle->SetPalette(55);
        hh->Draw("colz");
        hh->Draw("text same");

        corrCanv->Print("results/Correlation_histogram.C");
        corrCanv->Print("results/Correlation_histogram.pdf");
        corrCanv->Print("results/Correlation_histogram.eps");

    }



}   
