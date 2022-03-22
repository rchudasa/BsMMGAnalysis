
void MVATrainer::SetupTMVAOptions()
{
  std::cout<<"Invoking  SetupTMVAOptions() \n";
  Use["Cuts"]            = 0;
  Use["CutsD"]           = 0;
 
  Use["CutsPCA"]         = 0;
  Use["CutsGA"]          = 0;
  Use["CutsSA"]          = 0;
  
  // 1-dimensional likelihood ("naive Bayes estimator")
  Use["Likelihood"]      = 0;
  Use["LikelihoodD"]     = 0; // the "D" extension indicates decorrelated input variables (see option strings)
  Use["LikelihoodPCA"]   = 0; // the "PCA" extension indicates PCA-transformed input variables (see option strings)
  Use["LikelihoodKDE"]   = 0;
  Use["LikelihoodMIX"]   = 0;
  
  // Mutidimensional likelihood and Nearest-Neighbour methods
  Use["PDERS"]           = 0;
  Use["PDERSD"]          = 0;
  Use["PDERSPCA"]        = 0;
  Use["PDEFoam"]         = 0;
  Use["PDEFoamBoost"]    = 0; // uses generalised MVA method boosting
  Use["KNN"]             = 0; // k-nearest neighbour method
  
  // Linear Discriminant Analysis
  Use["LD"]              = 0; // Linear Discriminant identical to Fisher
  Use["Fisher"]          = 0;
  Use["FisherG"]         = 0;
  Use["BoostedFisher"]   = 0; // uses generalised MVA method boosting
  Use["HMatrix"]         = 0;
  
  // Function Discriminant analysis
  Use["FDA_GA"]          = 0; // minimisation of user-defined function using Genetics Algorithm
  Use["FDA_SA"]          = 0;
  Use["FDA_MC"]          = 0;
  Use["FDA_MT"]          = 0;
  Use["FDA_GAMT"]        = 0;
  Use["FDA_MCMT"]        = 0;
  
  // Neural Networks (all are feed-forward Multilayer Perceptrons)
  Use["MLP"]             = 0; // Recommended ANN
  Use["MLPBFGS"]         = 0; // Recommended ANN with optional training method
  Use["MLPBNN"]          = 0; // Recommended ANN with BFGS training method and bayesian regulator
  Use["CFMlpANN"]        = 0; // Depreciated ANN from ALEPH
  Use["TMlpANN"]         = 0; // ROOT's own ANN
  Use["DNN_GPU"]         = 0; // CUDA-accelerated DNN training.
  Use["DNN_CPU"]         = 0; // Multi-core accelerated DNN.
  
  // Support Vector Machine
  Use["SVM"]             = 0;
  
  // Boosted Decision Trees
  Use["BDT"]             = 0; // uses Adaptive Boost
  Use["BDTG"]            = 0; // uses Gradient Boost
  Use["BDTB"]            = 0; // uses Bagging
  Use["BDTD"]            = 0; // decorrelation + Adaptive Boost
  Use["BDTF"]            = 0; // allow usage of fisher discriminant for node splitting
  
  // Friedman's RuleFit method, ie, an optimised series of cuts ("rules")
  Use["RuleFit"]         = 0;

  // ---------------------------------------------------------------
  std::cout << std::endl;
  std::cout << "==> Start TMVAClassification" << std::endl;
    std::cout<<__LINE__<<" \n";
  
  // Select methods (don't look at this code - not of interest)
  if (mvaMethords != "") {
    for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;
    
    std::vector<TString> mlist = TMVA::gTools().SplitString( mvaMethords, ',' );
    for (UInt_t i=0; i<mlist.size(); i++) {
      std::string regMethod(mlist[i]);
      
      if (Use.find(regMethod) == Use.end()) {
	std::cout << "Method \"" << regMethod << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
	for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) std::cout << it->first << " ";
	std::cout << std::endl;
      return;
      }
      Use[regMethod] = 1;
    }
  }
  
  // --------------------------------------------------------------------------------------------------
  
  // Here the preparation phase begins
  
  // Read training and test data
  // (it is also possible to use ASCII format as input -> see TMVA Users Guide)
    std::cout<<__LINE__<<" \n";
  for(int i =0;i<mvaTrainVars.size();i++)
  {
        // (var, prettyTitle , unit ,type ) 
        std::cout<<"Adding "<<mvaTrainVars[i]<<" to "<<" data loader \n";
        dataloader->AddVariable(mvaTrainVars[i].c_str(),mvaTrainVars[i].c_str(),"",'F');
  }
  
  // You can add so-called "Spectator variables", which are not used in the MVA training,
  // but will appear in the final "TestTree" produced by TMVA. This TestTree will contain the
  // input variables, the response values of all trained MVAs, and the spectator variables
  for(int i =0;i<spectatorVars.size();i++)
  {
        // (var, prettyTitle , unit ,type ) 
        std::cout<<"Adding sectator "<<spectatorVars[i]<<" to "<<" data loader \n";
        dataloader->AddSpectator(spectatorVars[i].c_str(),spectatorVars[i].c_str(),"",'F');
  }
  
  
   signalFiles.clear();
   signalTrees.clear();


  // global event weights per tree (see below for setting event-wise weights)
  // You can add an arbitrary number of signal or background trees
    std::cout<<__LINE__<<" \n";
  Double_t x;
  for(int i=0 ; i < signalFnames.size() ; i++)
  {
    std::cout<<__LINE__<<" \n";
            auto aFile = TFile::Open( signalFnames[i].c_str() );
            signalFiles.push_back(aFile)     ;
            auto aTree =  (TTree * ) signalFiles[i]->Get(signalTreeNames[i].c_str());
            signalTrees.push_back(aTree);
            std::cout<<"Adding sig "<<signalTreeNames[i]<<" from "<<signalFnames[i]<<" with "<<signalTrees[i]->GetEntries()<<"\n";
            dataloader->AddSignalTree    ( signalTrees[i], signalWeight[i] );
            //for(int j =0;j<mvaTrainVars.size();j++)
            //{
            //    signalTrees[i]->SetBranchAddress(mvaTrainVars[j].c_str(),&x);
            //    std::cout<<" j =" <<j<<"\n";
            //}   
            //for(int j =0;j<spectatorVars.size();j++)
            //{
            //    signalTrees[i]->SetBranchAddress(spectatorVars[j].c_str(),&x);
            //    std::cout<<" j =" <<j<<"\n";
            //}   
  }
  for(int i=0 ; i < bkgFnames.size() ; i++)
  {
    std::cout<<__LINE__<<" \n";
            auto aFile = TFile::Open( bkgFnames[i].c_str() );
            bkgFiles.push_back(aFile)     ;
            auto aTree =  (TTree * ) bkgFiles[i]->Get(bkgTreeNames[i].c_str());
            bkgTrees.push_back(aTree);
            std::cout<<"Adding bkg "<<bkgTreeNames[i]<<" from "<<bkgFnames[i]<<" with "<<bkgTrees[i]->GetEntries()<<"\n";
            dataloader->AddBackgroundTree    ( bkgTrees[i], bkgWeight[i] );
            //for(int j =0;j<mvaTrainVars.size();j++)
            //{
            //    bkgTrees[i]->SetBranchAddress(mvaTrainVars[j].c_str(),&x);
            //    std::cout<<" j =" <<j<<"\n";
            //}  
            //for(int j =0;j<spectatorVars.size();j++)
            //{
            //    bkgTrees[i]->SetBranchAddress(spectatorVars[j].c_str(),&x);
            //    std::cout<<" j =" <<j<<"\n";
            //}  
    }


  // cuts can be like :  "VBFJet_mjj > 400 && leadVBF_pt > 40 && subleadVBF_pt > 30 "; // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
  // cuts can be like :  "VBFJet_mjj > 400 && leadVBF_pt > 40 && subleadVBF_pt > 30 "; // for example: TCut mycutb = "abs(var1)<0.5";
  dataloader->PrepareTrainingAndTestTree( sigCuts.c_str() , bkgCuts.c_str(), testTrainConfig.c_str() );
  std::cout<<sigCuts.c_str()<<" , "<<bkgCuts.c_str()<<" , "<<testTrainConfig<<"\n";
  
    std::cout<<__LINE__<<" \n";
  
  // ### Book MVA methods
  //
  // Please lookup the various method configuration options in the corresponding cxx files, eg:
  // src/MethoCuts.cxx, etc, or here: http://tmva.sourceforge.net/optionRef.html
  // it is possible to preset ranges in the option string in which the cut optimisation should be done:
  // "...:CutRangeMin[2]=-1:CutRangeMax[2]=1"...", where [2] is the third input variable
  
  // Cut optimisation

}

void MVATrainer::trainAndTest()
{
    std::cout<<__LINE__<<" \n";

  if (Use["Cuts"]) { 
     if(mvaMethordOptions["Cuts"]=="auto")
     {
            mvaMethordOptions["Cuts"]="!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart";
     }
    std::cout<<"Setting the mva option for "<<"CUTS"<<" as : "<<mvaMethordOptions["Cuts"]<<"\n";
    factory->BookMethod( dataloader, TMVA::Types::kCuts, "Cuts",mvaMethordOptions["Cuts"].c_str() ); 	}
  
  if (Use["CutsD"]) { 
    factory->BookMethod( dataloader, TMVA::Types::kCuts, "CutsD",
			 "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart:VarTransform=Decorrelate" );
	}
  
  if (Use["CutsPCA"]) { 
    factory->BookMethod( dataloader, TMVA::Types::kCuts, "CutsPCA",
			 "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart:VarTransform=PCA" );
	}
  
  if (Use["CutsGA"]) { 
    factory->BookMethod( dataloader, TMVA::Types::kCuts, "CutsGA",
			 "H:!V:FitMethod=GA:CutRangeMin[0]=-10:CutRangeMax[0]=10:VarProp[1]=FMax:EffSel:Steps=30:Cycles=3:PopSize=400:SC_steps=10:SC_rate=5:SC_factor=0.95" );
	}
  
  if (Use["CutsSA"]) { 
    factory->BookMethod( dataloader, TMVA::Types::kCuts, "CutsSA",
			 "!H:!V:FitMethod=SA:EffSel:MaxCalls=150000:KernelTemp=IncAdaptive:InitialTemp=1e+6:MinTemp=1e-6:Eps=1e-10:UseDefaultScale" );
	}
  
  // Likelihood ("naive Bayes estimator")
  if (Use["Likelihood"]) { 
    factory->BookMethod( dataloader, TMVA::Types::kLikelihood, "Likelihood",
			 "H:!V:TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmoothBkg[1]=10:NSmooth=1:NAvEvtPerBin=50" );
	}
  
  // Decorrelated likelihood
  if (Use["LikelihoodD"]) { 
    factory->BookMethod( dataloader, TMVA::Types::kLikelihood, "LikelihoodD",
			 "!H:!V:TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmooth=5:NAvEvtPerBin=50:VarTransform=Decorrelate" );
	}
  
  // PCA-transformed likelihood
  if (Use["LikelihoodPCA"]) { 
    factory->BookMethod( dataloader, TMVA::Types::kLikelihood, "LikelihoodPCA",
			 "!H:!V:!TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmooth=5:NAvEvtPerBin=50:VarTransform=PCA" );
	}
  
  // Use a kernel density estimator to approximate the PDFs
  if (Use["LikelihoodKDE"]) { 
    factory->BookMethod( dataloader, TMVA::Types::kLikelihood, "LikelihoodKDE",
			 "!H:!V:!TransformOutput:PDFInterpol=KDE:KDEtype=Gauss:KDEiter=Adaptive:KDEFineFactor=0.3:KDEborder=None:NAvEvtPerBin=50" );
	}
  
  // Use a variable-dependent mix of splines and kernel density estimator
  if (Use["LikelihoodMIX"]) { 
    factory->BookMethod( dataloader, TMVA::Types::kLikelihood, "LikelihoodMIX",
			 "!H:!V:!TransformOutput:PDFInterpolSig[0]=KDE:PDFInterpolBkg[0]=KDE:PDFInterpolSig[1]=KDE:PDFInterpolBkg[1]=KDE:PDFInterpolSig[2]=Spline2:PDFInterpolBkg[2]=Spline2:PDFInterpolSig[3]=Spline2:PDFInterpolBkg[3]=Spline2:KDEtype=Gauss:KDEiter=Nonadaptive:KDEborder=None:NAvEvtPerBin=50" );
	}
  
  // Test the multi-dimensional probability density estimator
  // here are the options strings for the MinMax and RMS methods, respectively:
  //
  //      "!H:!V:VolumeRangeMode=MinMax:DeltaFrac=0.2:KernelEstimator=Gauss:GaussSigma=0.3" );
  //      "!H:!V:VolumeRangeMode=RMS:DeltaFrac=3:KernelEstimator=Gauss:GaussSigma=0.3" );
  if (Use["PDERS"]) { 
    factory->BookMethod( dataloader, TMVA::Types::kPDERS, "PDERS",
			 "!H:!V:NormTree=T:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600" );
	}
  
  if (Use["PDERSD"]) { 
    factory->BookMethod( dataloader, TMVA::Types::kPDERS, "PDERSD",
			 "!H:!V:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600:VarTransform=Decorrelate" );
	}
  
  if (Use["PDERSPCA"]) { 
    factory->BookMethod( dataloader, TMVA::Types::kPDERS, "PDERSPCA",
			 "!H:!V:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600:VarTransform=PCA" );
	}
  
  // Multi-dimensional likelihood estimator using self-adapting phase-space binning
  if (Use["PDEFoam"]) { 
    factory->BookMethod( dataloader, TMVA::Types::kPDEFoam, "PDEFoam",
			 "!H:!V:SigBgSeparate=F:TailCut=0.001:VolFrac=0.0666:nActiveCells=500:nSampl=2000:nBin=5:Nmin=100:Kernel=None:Compress=T" );
	}
  
  if (Use["PDEFoamBoost"]) { 
    factory->BookMethod( dataloader, TMVA::Types::kPDEFoam, "PDEFoamBoost",
			 "!H:!V:Boost_Num=30:Boost_Transform=linear:SigBgSeparate=F:MaxDepth=4:UseYesNoCell=T:DTLogic=MisClassificationError:FillFoamWithOrigWeights=F:TailCut=0:nActiveCells=500:nBin=20:Nmin=400:Kernel=None:Compress=T" );
	}
  
  // K-Nearest Neighbour classifier (KNN)
  if (Use["KNN"]) { 
    factory->BookMethod( dataloader, TMVA::Types::kKNN, "KNN",
			 "H:nkNN=20:ScaleFrac=0.8:SigmaFact=1.0:Kernel=Gaus:UseKernel=F:UseWeight=T:!Trim" );
	}
  
  // H-Matrix (chi2-squared) method
  if (Use["HMatrix"]) { 
    factory->BookMethod( dataloader, TMVA::Types::kHMatrix, "HMatrix", "!H:!V:VarTransform=None" );
	}
  
  // Linear discriminant (same as Fisher discriminant)
  if (Use["LD"]) { 
    factory->BookMethod( dataloader, TMVA::Types::kLD, "LD", "H:!V:VarTransform=None:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10" );
	}
  
  // Fisher discriminant (same as LD)
  if (Use["Fisher"]) { 
    factory->BookMethod( dataloader, TMVA::Types::kFisher, "Fisher", "H:!V:Fisher:VarTransform=None:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10" );
	}
  
  // Fisher with Gauss-transformed input variables
  if (Use["FisherG"]) { 
    factory->BookMethod( dataloader, TMVA::Types::kFisher, "FisherG", "H:!V:VarTransform=Gauss" );
	}
  
  // Composite classifier: ensemble (tree) of boosted Fisher classifiers
  if (Use["BoostedFisher"]) { 
    factory->BookMethod( dataloader, TMVA::Types::kFisher, "BoostedFisher",
			 "H:!V:Boost_Num=20:Boost_Transform=log:Boost_Type=AdaBoost:Boost_AdaBoostBeta=0.2:!Boost_DetailedMonitoring" );
	}
  
  // Function discrimination analysis (FDA) -- test of various fitters - the recommended one is Minuit (or GA or SA)
  if (Use["FDA_MC"]) { 
    factory->BookMethod( dataloader, TMVA::Types::kFDA, "FDA_MC",
			 "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MC:SampleSize=100000:Sigma=0.1" );
	}
  
  if (Use["FDA_GA"]) {  // can also use Simulated Annealing (SA) algorithm (see Cuts_SA options])
    factory->BookMethod( dataloader, TMVA::Types::kFDA, "FDA_GA",
			 "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=GA:PopSize=100:Cycles=2:Steps=5:Trim=True:SaveBestGen=1" );
	}
  
  if (Use["FDA_SA"]) {  // can also use Simulated Annealing (SA) algorithm (see Cuts_SA options])
    factory->BookMethod( dataloader, TMVA::Types::kFDA, "FDA_SA",
			 "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=SA:MaxCalls=15000:KernelTemp=IncAdaptive:InitialTemp=1e+6:MinTemp=1e-6:Eps=1e-10:UseDefaultScale" );
	}
  
  if (Use["FDA_MT"]) { 
    factory->BookMethod( dataloader, TMVA::Types::kFDA, "FDA_MT",
			 "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=2:UseImprove:UseMinos:SetBatch" );
	}
  
  if (Use["FDA_GAMT"]) { 
    factory->BookMethod( dataloader, TMVA::Types::kFDA, "FDA_GAMT",
			 "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=GA:Converger=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=0:!UseImprove:!UseMinos:SetBatch:Cycles=1:PopSize=5:Steps=5:Trim" );
	}
  
  if (Use["FDA_MCMT"]) { 
    factory->BookMethod( dataloader, TMVA::Types::kFDA, "FDA_MCMT",
			 "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MC:Converger=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=0:!UseImprove:!UseMinos:SetBatch:SampleSize=20" );
	}
  
  // TMVA ANN: MLP (recommended ANN) -- all ANNs in TMVA are Multilayer Perceptrons
  if (Use["MLP"]) { 
    factory->BookMethod( dataloader, TMVA::Types::kMLP, "MLP", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=100:HiddenLayers=N+5:TestRate=5:!UseRegulator" );
	}
  
  if (Use["MLPBFGS"]) { 
    factory->BookMethod( dataloader, TMVA::Types::kMLP, "MLPBFGS", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:!UseRegulator" );
	}
  
  if (Use["MLPBNN"]) { 
    factory->BookMethod( dataloader, TMVA::Types::kMLP, "MLPBNN", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=60:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:UseRegulator" );
	} // BFGS training with bayesian regulators
  
  
  // Multi-architecture DNN implementation.


    if (Use["DNN_CPU"] or Use["DNN_GPU"]) {  

     TString layoutString ("Layout=TANH|128,TANH|128,TANH|128,LINEAR");
     if(Use["DNN_CPU"]) {
        if(mvaMethordOptions["DNN_CPU"]=="auto")
        {
               mvaMethordOptions["DNN_CPU"]="Layout=TANH|128,TANH|128,TANH|128,LINEAR";
        }  // Gradient Boost
        std::cout<<"Setting the mva option for "<<"DNN_CPU"<<" as : "<<mvaMethordOptions["DNN_CPU"]<<"\n";
        layoutString=mvaMethordOptions["DNN_CPU"].c_str();
     }

     if(Use["DNN_GPU"]) {
        if(mvaMethordOptions["DNN_GPU"]=="auto")
        {
               mvaMethordOptions["DNN_GPU"]="Layout=TANH|128,TANH|128,TANH|128,LINEAR";
        }  // Gradient Boost
        std::cout<<"Setting the mva option for "<<"DNN_GPU"<<" as : "<<mvaMethordOptions["DNN_GPU"]<<"\n";
        layoutString=mvaMethordOptions["DNN_GPU"].c_str();
     }


     // General layout.
     // Training strategies.
     TString training0("LearningRate=1e-1,Momentum=0.9,Repetitions=1,"
		       "ConvergenceSteps=20,BatchSize=256,TestRepetitions=10,"
		       "WeightDecay=1e-4,Regularization=L2,"
		       "DropConfig=0.0+0.5+0.5+0.5, Multithreading=True");
     TString training1("LearningRate=1e-2,Momentum=0.9,Repetitions=1,"
		       "ConvergenceSteps=20,BatchSize=256,TestRepetitions=10,"
		       "WeightDecay=1e-4,Regularization=L2,"
		       "DropConfig=0.0+0.0+0.0+0.0, Multithreading=True");
     TString training2("LearningRate=1e-3,Momentum=0.0,Repetitions=1,"
		       "ConvergenceSteps=20,BatchSize=256,TestRepetitions=10,"
		       "WeightDecay=1e-4,Regularization=L2,"
		       "DropConfig=0.0+0.0+0.0+0.0, Multithreading=True");
     TString trainingStrategyString ("TrainingStrategy=");
     trainingStrategyString += training0 + "|" + training1 + "|" + training2;
     
     // General Options.
     TString dnnOptions ("!H:V:ErrorStrategy=CROSSENTROPY:VarTransform=N:"
			 "WeightInitialization=XAVIERUNIFORM");
     dnnOptions.Append (":"); dnnOptions.Append (layoutString);
     dnnOptions.Append (":"); dnnOptions.Append (trainingStrategyString);
     
     // Cuda implementation.
     if (Use["DNN_GPU"]) {  
       TString gpuOptions = dnnOptions + ":Architecture=GPU";
       factory->BookMethod(dataloader, TMVA::Types::kDNN, "DNN_GPU", gpuOptions);
     }
     // Multi-core CPU implementation.
     if (Use["DNN_CPU"]) {  
       TString cpuOptions = dnnOptions + ":Architecture=CPU";
       factory->BookMethod(dataloader, TMVA::Types::kDNN, "DNN_CPU", cpuOptions);
     }
   }
   
   // CF(Clermont-Ferrand)ANN
   if (Use["CFMlpANN"]) { 
     factory->BookMethod( dataloader, TMVA::Types::kCFMlpANN, "CFMlpANN", "!H:!V:NCycles=200:HiddenLayers=N+1,N"  );
	} // n_cycles:#nodes:#nodes:...
   
   // Tmlp(Root)ANN
   if (Use["TMlpANN"]) { 
     factory->BookMethod( dataloader, TMVA::Types::kTMlpANN, "TMlpANN", "!H:!V:NCycles=200:HiddenLayers=N+1,N:LearningMethod=BFGS:ValidationFraction=0.3"  );
	} // n_cycles:#nodes:#nodes:...
   
   // Support Vector Machine
   if (Use["SVM"]) { 
     factory->BookMethod( dataloader, TMVA::Types::kSVM, "SVM", "Gamma=0.25:Tol=0.001:VarTransform=Norm" );
	}
   
   // Boosteda Decision Trees
   if (Use["BDTG"]) { 
     if(mvaMethordOptions["BDTG"]=="auto")
     {
            mvaMethordOptions["BDTG"]="!H:!V:NTrees=1000:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=3";

     }  // Gradient Boost
     std::cout<<"Setting the mva option for "<<"BDTG"<<" as : "<<mvaMethordOptions["BDTG"]<<"\n";
     factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDTG", mvaMethordOptions["BDTG"].c_str() );
	}
   
   if (Use["BDT"]) {  
    // Adaptive Boost
     if(mvaMethordOptions["BDT"]=="auto")
     {
       mvaMethordOptions["BDT"]="!H:V=True:NTrees=500:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20" ;
     }     
     std::cout<<"Setting the mva option for "<<"BDT"<<" as : "<<mvaMethordOptions["BDT"]<<"\n";
     factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDT",mvaMethordOptions["BDT"]	);
	}
    
   if (Use["BDTB"]) {  // Bagging
     if(mvaMethordOptions["BDTB"]=="auto")
     {
       mvaMethordOptions["BDTB"]="!H:!V:NTrees=400:BoostType=Bagging:SeparationType=GiniIndex:nCuts=20"  ;
     }     
     std::cout<<"Setting the mva option for "<<"CUTS"<<" as : "<<mvaMethordOptions["Cuts"]<<"\n";
     factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDTB", mvaMethordOptions["BDTB"].c_str() );
	}
   
   if (Use["BDTD"]) {  // Decorrelation + Adaptive Boost
     if(mvaMethordOptions["BDTD"]=="auto")
     {
       mvaMethordOptions["BDTD"] ="!H:!V:NTrees=400:MinNodeSize=5%:MaxDepth=3:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:VarTransform=Decorrelate"  ;     
     }      
     std::cout<<"Setting the mva option for "<<"BDTD"<<" as : "<<mvaMethordOptions["BDTD"]<<"\n";
     factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDTD", mvaMethordOptions["BDTD"].c_str());
	}
   
   if (Use["BDTF"]) {   // Allow Using Fisher discriminant in node splitting for (strong) linearly correlated variables
     
     if(mvaMethordOptions["BDTF"]=="auto")
     {
       mvaMethordOptions["BDTF"]="!H:!V:NTrees=50:MinNodeSize=2.5%:UseFisherCuts:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20";     
     }     
     std::cout<<"Setting the mva option for "<<"BDTF"<<" as : "<<mvaMethordOptions["BDTF"]<<"\n";
     factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDTF", mvaMethordOptions["BDTF"] );
	}
   
   // RuleFit -- TMVA implementation of Friedman's method
   if (Use["RuleFit"]) { 
     factory->BookMethod( dataloader, TMVA::Types::kRuleFit, "RuleFit",
			  "H:!V:RuleFitModule=RFTMVA:Model=ModRuleLinear:MinImp=0.001:RuleMinDist=0.001:NTrees=20:fEventsMin=0.01:fEventsMax=0.5:GDTau=-1.0:GDTauPrec=0.01:GDStep=0.01:GDNSteps=10000:GDErrScale=1.02" );
	}
   
   // For an example of the category classifier usage, see: TMVAClassificationCategory
   //
   // --------------------------------------------------------------------------------------------------
   //  Now you can optimize the setting (configuration) of the MVAs using the set of training events
   // STILL EXPERIMENTAL and only implemented for BDT's !
   //
   //     factory->OptimizeAllMethods("SigEffAt001","Scan");
   //     factory->OptimizeAllMethods("ROCIntegral","FitGA");
   //
   // --------------------------------------------------------------------------------------------------
   
   // Now you can tell the factory to train, test, and evaluate the MVAs
   //
   // Train MVAs using the set of training events
   std::cout<<__LINE__<<"\n";
   factory->TrainAllMethods();
   
   // Evaluate all MVAs using the set of test events
   factory->TestAllMethods();
   
   // Evaluate and compare performance of all configured MVAs
   factory->EvaluateAllMethods();
   
   // --------------------------------------------------------------
   
   
   std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
   std::cout << "==> TMVAClassification is done!" << std::endl;
   
   // Launch the GUI for the root macros
   //if (!gROOT->IsBatch()) TMVA::TMVAGui( outfileName );
   //   Double_t 	area= TMVA::ROCCalc::GetROCIntegral();
   //   std::cout << area << endl; 	
   
}

