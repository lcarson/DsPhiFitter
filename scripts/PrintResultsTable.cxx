void PrintResultsTable(std::string filelocation = "/home/hadavizadeh/DST/B2DsPhi/Systematics/2017-03-14_Nominal_fit_Merged_DsPhi/Bu2DsPhi_Fitter_copy/systematicsDir/vary_lastRooFitResult_summed_splitHel_s21_s21r1_s24_s26_Ds2PhiPi_Ds2KKPi_Ds2PiPiPi_Ds2KPiPi_SEED_700_16243.root"){
    
    filelocation = "/home/hadavizadeh/DST/B2DsPhi/Systematics/2017-03-26_Nominal_unmerged_DsPhi/Bu2DsPhi_Fitter_copy/systematicsDir/vary_lastRooFitResult_summed_splitHel_s21_s21r1_s24_s26_Ds2PhiPi_Ds2KKPi_Ds2PiPiPi_Ds2KPiPi_SEED_300_13439.root";

    std::map<std::string,std::string> Description;
    Description["Branching_fraction"]           = "Branching fraction ($\\times10^{-7}$)                             ";
    Description["Dsa1_to_DsstPhi_fraction"]     = "Fraction of $D_{s}^{*} \\phi$ in ($D_{s}^{(*)}KK^{*0}$ + $D_{s}^{*} \\phi$)";
    Description["Fraction_DsstD0_Helbin1"]      = "Fraction of $D_{s}^{*}D^{0}$ in ($D_{s}^{*}D^{0}$+$D_{s}D^{*0}$) in 1st helicity bin  ";
    Description["Fraction_DsstD0_Helbin2"]      = "Fraction of $D_{s}^{*}D^{0}$ in ($D_{s}^{*}D^{0}$+$D_{s}D^{*0}$) in 2nd helicity bin  ";
    Description["Ratio_DD_to_Dsa1"]             = "Ratio of $DD^{'}$ backgrounds to ($D_{s}^{*} \\phi$ + $D_{s}^{(*)}KK^{*0}$)";
    Description["Ratio_DsstDst0_to_Lowmass"]    = "Ratio of $D_{s}^{*}D^{*0}$ to ($D_{s}^{*}D^{0}$ + $D_{s}D^{*0}$)";
    Description["global_comb_slope"]            = "Combinatorial slope                                             ";
    Description["frac_Dsa1_DsstKKst"]           = "Fraction of $D_{s}KK^{*0}$ in ($D_{s}^{*}KK^{*0}$+$D_{s}KK^{*0}$)";
    Description["global_csi"]                   = "Relative heights double peaked part reco shapes                 ";
    Description["global_mean"]                  = "Mean B mass (\\mev)                                              ";
    Description["global_shift"]                 = "Part reco background mass offset (\\mev)                         ";
    Description["low_mass_total_DsD0_Ds2KKPi_both_both_DsBDTbin1_PhiBDTbin1_both_both"]    = "Yield of ($D_{s}^{*}D^{0}$ + $D_{s}D^{*0}$) in $D_{s} \\to KK\\pi$";
    Description["low_mass_total_DsD0_Ds2KPiPi_both_both_DsBDTbin1_PhiBDTbin1_both_both"]   = "Yield of ($D_{s}^{*}D^{0}$ + $D_{s}D^{*0}$) in $D_{s} \\to K\\pi\\pi$";
    Description["low_mass_total_DsD0_Ds2PhiPi_both_both_DsBDTbin1_PhiBDTbin1_both_both"]   = "Yield of ($D_{s}^{*}D^{0}$ + $D_{s}D^{*0}$) in $D_{s} \\to \\phi\\pi$";
    Description["low_mass_total_DsD0_Ds2PiPiPi_both_both_DsBDTbin1_PhiBDTbin1_both_both"]  = "Yield of ($D_{s}^{*}D^{0}$ + $D_{s}D^{*0}$) in $D_{s} \\to \\pi\\pi\\pi$";
    Description["low_mass_total_DsPhi_Ds2KKPi_both_both_DsBDTbin1_PhiBDTbin1_both_both"]   = "Yield of ($D_{s}^{*} \\phi$ + $D_{s}^{(*)}KK^{*0}$) in $D_{s} \\to KK\\pi$";
    Description["low_mass_total_DsPhi_Ds2KPiPi_both_both_DsBDTbin1_PhiBDTbin1_both_both"]  = "Yield of ($D_{s}^{*} \\phi$ + $D_{s}^{(*)}KK^{*0}$) in $D_{s} \\to K\\pi\\pi$";
    Description["low_mass_total_DsPhi_Ds2PhiPi_both_both_DsBDTbin1_PhiBDTbin1_both_both"]  = "Yield of ($D_{s}^{*} \\phi$ + $D_{s}^{(*)}KK^{*0}$) in $D_{s} \\to \\phi\\pi$";
    Description["low_mass_total_DsPhi_Ds2PiPiPi_both_both_DsBDTbin1_PhiBDTbin1_both_both"] = "Yield of ($D_{s}^{*} \\phi$ + $D_{s}^{(*)}KK^{*0}$) in $D_{s} \\to K\\pi\\pi$";
    Description["sigma_DsD0_Ds2KKPi"]           = "$\\sigma_{1}$ for $D_{s} \\to KK\\pi$  (\\mev)                          ";
    Description["sigma_DsD0_Ds2KPiPi"]          = "$\\sigma_{1}$ for $D_{s} \\to K\\pi\\pi$ (\\mev)                        ";
    Description["sigma_DsD0_Ds2PhiPi"]          = "$\\sigma_{1}$ for $D_{s} \\to \\phi\\pi$  (\\mev)                       ";
    Description["sigma_DsD0_Ds2PiPiPi"]         = "$\\sigma_{1}$ for $D_{s} \\to \\pi\\pi\\pi$  (\\mev)                    ";
    Description["splitHel_DsD0_PR_peak_fraction"] = "Fraction of normalisation part reco in 1st helicity bin";
    Description["splitHel_DsD0_peak_fraction"]    = "Fraction of normalisation peak in 1st helicity bin";
    
    Description["yield_comb_DsD0_Ds2KKPi_both_Helbin1_DsBDTbin1_PhiBDTbin1_both_both"]        = "Comb. Yield $D_{s}D^{0}$ $D_{s} \\to KK\\pi$ in 1st helicity bin";
    Description["yield_comb_DsD0_Ds2KKPi_both_Helbin2_DsBDTbin1_PhiBDTbin1_both_both"]        = "Comb. Yield $D_{s}D^{0}$ $D_{s} \\to KK\\pi$ in 2nd helicity bin";
    Description["yield_comb_DsD0_Ds2KPiPi_both_Helbin1_DsBDTbin1_PhiBDTbin1_both_both"]       = "Comb. Yield $D_{s}D^{0}$ $D_{s} \\to K\\pi\\pi$ in 1st helicity bin";
    Description["yield_comb_DsD0_Ds2KPiPi_both_Helbin2_DsBDTbin1_PhiBDTbin1_both_both"]       = "Comb. Yield $D_{s}D^{0}$ $D_{s} \\to K\\pi\\pi$ in 2nd helicity bin";
    Description["yield_comb_DsD0_Ds2PhiPi_both_Helbin1_DsBDTbin1_PhiBDTbin1_both_both"]       = "Comb. Yield $D_{s}D^{0}$ $D_{s} \\to \\phi\\pi$ in 1st helicity bin";
    Description["yield_comb_DsD0_Ds2PhiPi_both_Helbin2_DsBDTbin1_PhiBDTbin1_both_both"]       = "Comb. Yield $D_{s}D^{0}$ $D_{s} \\to \\phi\\pi$ in 2nd helicity bin";
    Description["yield_comb_DsD0_Ds2PiPiPi_both_Helbin1_DsBDTbin1_PhiBDTbin1_both_both"]      = "Comb. Yield $D_{s}D^{0}$ $D_{s} \\to \\pi\\pi\\pi$ in 1st helicity bin";
    Description["yield_comb_DsD0_Ds2PiPiPi_both_Helbin2_DsBDTbin1_PhiBDTbin1_both_both"]      = "Comb. Yield $D_{s}D^{0}$ $D_{s} \\to \\pi\\pi\\pi$ in 2nd helicity bin";
    
    Description["yield_comb_DsPhiSide_Ds2KKPi_both_Helbin1_DsBDTbin1_PhiBDTbin1_both_both"]   = "Comb. Yield $D_{s}\\phi$ sideband $D_{s} \\to KK\\pi$ in 1st helicity bin";
    Description["yield_comb_DsPhiSide_Ds2KKPi_both_Helbin2_DsBDTbin1_PhiBDTbin1_both_both"]   = "Comb. Yield $D_{s}\\phi$ sideband $D_{s} \\to KK\\pi$ in 2nd helicity bin";
    Description["yield_comb_DsPhiSide_Ds2KPiPi_both_Helbin1_DsBDTbin1_PhiBDTbin1_both_both"]  = "Comb. Yield $D_{s}\\phi$ sideband $D_{s} \\to K\\pi\\pi$ in 1st helicity bin";
    Description["yield_comb_DsPhiSide_Ds2KPiPi_both_Helbin2_DsBDTbin1_PhiBDTbin1_both_both"]  = "Comb. Yield $D_{s}\\phi$ sideband $D_{s} \\to K\\pi\\pi$ in 2nd helicity bin";
    Description["yield_comb_DsPhiSide_Ds2PhiPi_both_Helbin1_DsBDTbin1_PhiBDTbin1_both_both"]  = "Comb. Yield $D_{s}\\phi$ sideband $D_{s} \\to \\phi\\pi$ in 1st helicity bin";
    Description["yield_comb_DsPhiSide_Ds2PhiPi_both_Helbin2_DsBDTbin1_PhiBDTbin1_both_both"]  = "Comb. Yield $D_{s}\\phi$ sideband $D_{s} \\to \\phi\\pi$ in 2nd helicity bin";
    Description["yield_comb_DsPhiSide_Ds2PiPiPi_both_Helbin1_DsBDTbin1_PhiBDTbin1_both_both"] = "Comb. Yield $D_{s}\\phi$ sideband $D_{s} \\to \\pi\\pi\\pi$ in 1st helicity bin";
    Description["yield_comb_DsPhiSide_Ds2PiPiPi_both_Helbin2_DsBDTbin1_PhiBDTbin1_both_both"] = "Comb. Yield $D_{s}\\phi$ sideband $D_{s} \\to \\pi\\pi\\pi$ in 2nd helicity bin";
    
    Description["yield_comb_DsPhi_Ds2KKPi_both_Helbin1_DsBDTbin1_PhiBDTbin1_both_both"]       = "Comb. Yield $D_{s}\\phi$ $D_{s} \\to KK\\pi$ in 1st helicity bin";
    Description["yield_comb_DsPhi_Ds2KKPi_both_Helbin2_DsBDTbin1_PhiBDTbin1_both_both"]       = "Comb. Yield $D_{s}\\phi$ $D_{s} \\to KK\\pi$ in 2nd helicity bin";
    Description["yield_comb_DsPhi_Ds2KPiPi_both_Helbin1_DsBDTbin1_PhiBDTbin1_both_both"]      = "Comb. Yield $D_{s}\\phi$ $D_{s} \\to K\\pi\\pi$ in 1st helicity bin";
    Description["yield_comb_DsPhi_Ds2KPiPi_both_Helbin2_DsBDTbin1_PhiBDTbin1_both_both"]      = "Comb. Yield $D_{s}\\phi$ $D_{s} \\to K\\pi\\pi$ in 2nd helicity bin";
    Description["yield_comb_DsPhi_Ds2PhiPi_both_Helbin1_DsBDTbin1_PhiBDTbin1_both_both"]      = "Comb. Yield $D_{s}\\phi$ $D_{s} \\to \\phi\\pi$ in 1st helicity bin";
    Description["yield_comb_DsPhi_Ds2PhiPi_both_Helbin2_DsBDTbin1_PhiBDTbin1_both_both"]      = "Comb. Yield $D_{s}\\phi$ $D_{s} \\to \\phi\\pi$ in 2nd helicity bin";
    Description["yield_comb_DsPhi_Ds2PiPiPi_both_Helbin1_DsBDTbin1_PhiBDTbin1_both_both"]     = "Comb. Yield $D_{s}\\phi$ $D_{s} \\to \\pi\\pi\\pi$ in 1st helicity bin";
    Description["yield_comb_DsPhi_Ds2PiPiPi_both_Helbin2_DsBDTbin1_PhiBDTbin1_both_both"]     = "Comb. Yield $D_{s}\\phi$ $D_{s} \\to \\pi\\pi\\pi$ in 2nd helicity bin";
    
    Description["yield_peak_DsD0_Ds2KKPi_both_both_DsBDTbin1_PhiBDTbin1_both_both"]           = "Total Normalisation Yield $D_{s}D^{0}$ $D_{s} \\to KK\\pi$ ";
    Description["yield_peak_DsD0_Ds2KPiPi_both_both_DsBDTbin1_PhiBDTbin1_both_both"]          = "Total Normalisation Yield $D_{s}D^{0}$ $D_{s} \\to K\\pi\\pi$ ";
    Description["yield_peak_DsD0_Ds2PhiPi_both_both_DsBDTbin1_PhiBDTbin1_both_both"]          = "Total Normalisation Yield $D_{s}D^{0}$ $D_{s} \\to \\phi\\pi$ ";
    Description["yield_peak_DsD0_Ds2PiPiPi_both_both_DsBDTbin1_PhiBDTbin1_both_both"]         = "Total Normalisation Yield $D_{s}D^{0}$ $D_{s} \\to \\pi\\pi\\pi$ ";


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

    std::cout << std::endl;
    std::cout << "\\begin{landscape}"<< std::endl;
    std::cout << "\\centering"                       << std::endl;
    std::cout << "\\begin{longtable}{|l|c|l|l|l|}"   << std::endl;
    std::cout << "\\hline"                           << std::endl;
    std::cout << "Variable & Value & Symmetric Error & Asym High Error & Asym Low Error \\\\"              << std::endl;
    std::cout << "\\hline"                           << std::endl;
    
    RooArgList arglist = result->floatParsFinal();

    RooRealVar* temp_var=0;
    int nv=0;
    int n_canvas=0;
    TIterator* it3 = arglist.createIterator();
    while((temp_var = (RooRealVar*)it3->Next())) {
      if(temp_var->isConstant() or temp_var->InheritsFrom("RooAbsCategory") ) continue;
      nv++;
      std::string varname = temp_var->GetName();
      std::cout << Description[varname] << "& "<< (varname=="Branching_fraction"?"xxx":Form("%f",temp_var->getVal())) <<" & $\\pm$ " << temp_var->getError()<< " & +" << temp_var->getAsymErrorHi()<< " & " << temp_var->getAsymErrorLo() <<" \\\\ " << std::endl;
                    
      //std::cout << "Name: " <<  temp_var->GetName() << "\t Description: "<< Description[temp_var->GetName()] << "\t"<<temp_var->getVal() << " +/- "<<temp_var->getError()<<std::endl;
    }


    std::cout << "\\hline"<< std::endl;
    std::cout << "\\caption{result} " << std::endl;
    std::cout << "\\label{tab:result} "<< std::endl;
    std::cout << "\\end{longtable}"<< std::endl;
    std::cout << "\\end{landscape}"<< std::endl;
    std::cout << std::endl;
}   
