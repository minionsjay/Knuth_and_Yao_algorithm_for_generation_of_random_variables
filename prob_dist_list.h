void which_prob_dist(int case_nb, struct prob_dist_data & T)
{
  /*
    0: simple 3 mass describe in the PDF

    1: discrete gaussian (this is not a discretization of the gaussian...)

    2: binomial

    3: a kind distribution arising from quantum mechanics 

    4: truncation of a zeta-dirichlet
    
    5: discretization of a singular distribution
  */
  if(case_nb==0)
    {
    
      T.nb_atoms = 3;

      NTL::RR PI = NTL::ComputePi_RR();
      NTL::RR EXminus1 = NTL::exp(NTL::RR(-1.0));
      
      T.pmv.push_back(NTL::RR(1.0)/PI);
      
      T.pmv.push_back(EXminus1);
      
      T.pmv.push_back(NTL::RR(1.0)-T.pmv[0]-T.pmv[1]);

    }
  else if(case_nb==1)
    {
      NTL::RR alpha = NTL::exp(NTL::RR(-1.0));
      NTL::RR beta  = NTL::RR(1000.0);
      
      NTL::ZZ zleft = NTL::conv<NTL::ZZ>(NTL::floor(alpha - NTL::RR(6.0)*beta));
      NTL::ZZ zright= NTL::conv<NTL::ZZ>(NTL::ceil(alpha + NTL::RR(6.0)*beta));
      
      NTL::RR norm_const = NTL::RR(0.0);
      T.nb_atoms = NTL::conv<long>(zright)-NTL::conv<long>(zleft)+1;
   
      NTL::ZZ z = zleft;
      
      do
	{
	  NTL::RR x = NTL::power( (NTL::conv<NTL::RR>(z)-alpha)/beta , 2.0);
	  NTL::RR y = NTL::exp(-x);
	  	  
	  norm_const = norm_const + y;
	  
	  T.pmv.push_back(y);
	  //std::cout << z << " --> " << y << "\n";
	  z+=NTL::ZZ(1);
       	}
      while(z<=zright);

      for(unsigned long i = 0 ; i < T.pmv.size() ; i++)
	{
	  T.pmv[i]=T.pmv[i]/norm_const;
	}
      
    }
  else if(case_nb==2)
    {
      long NB_TRIALS = 2000;
      NTL::RR SUCCESS_PROB = NTL::RR(0.15);

      long n = NB_TRIALS;
      NTL::RR p =SUCCESS_PROB;
      
      NTL::ZZ ** PT = new NTL::ZZ * [n+1];
      for(long i = 0 ; i < n+1; i++)
	{
	  PT[i] = new NTL::ZZ [i+1];
	}
      PT[0][0] = NTL::ZZ(1);
      PT[1][0] = NTL::ZZ(1);
      PT[1][1] = NTL::ZZ(1);
      
      for(long i = 2 ; i < n+1; i++)
	{
	  PT[i][0]=NTL::ZZ(1);
	  PT[i][i]=NTL::ZZ(1);
	  for(long j = 1; j<i ;j++)
	    {
	      PT[i][j] = PT[i-1][j-1]+PT[i-1][j];
	      
	    }

	}
      T.nb_atoms=n+1;
            
      for(long i=0;i<n+1;i++)
	{
	  T.pmv.push_back( NTL::conv<NTL::RR>(PT[n][i])*NTL::power(p,i)*NTL::power(1-p,n-i) );
	  delete [] PT[i];
	}
      delete [] PT;
    }
  else if(case_nb==3)
    {
      long n = 10;//dimension of Hilbert space = 2^n
      
      NTL::RR *theta = new NTL::RR [n];
      NTL::RR *phi   = new NTL::RR [n];

      T.nb_atoms=(1L<<n);
      
      NTL::RR PI = NTL::ComputePi_RR();
      //NTL::RR PI = NTL::conv<NTL::RR>("3.1415926535897932384626433832795029");
      for(long i = 0; i<n ; i++)
	{
	  
	  double f = 1.0/(2.0*n);
	  double g = 2*f;
	  NTL::RR sigm1 = NTL::conv<NTL::RR>(NTL::Jacobi(NTL::ZZ(i),NTL::ZZ((1L<<((n/2)-1)))));
	  for(long d = 0;d<n;d++)
	    {
	      if(i%n == d)
		{
		  theta[i] = sigm1*PI*NTL::RR(d*g+f);
		}
	    }
	
	  NTL::RR sigp1 = NTL::conv<NTL::RR>(NTL::Jacobi(NTL::ZZ(i),NTL::ZZ((1L<<((n/2)+1)))));
	  
	  for(long d = 0;d<n;d++)
	    {
	      if(i%n == d)
		{
		  phi[i] = sigp1*NTL::RR(d*g+f);
		}
	    }
	  
	}
      
      NTL::RR thetasum = NTL::RR(0.0);
      for(long i = 0;i<n;i++)
	{
	  thetasum+=theta[i];
	}
      thetasum = NTL::RR(0.5)*thetasum;
      std::cout << "theta = " << thetasum/PI << "\n";
      //thetasum = NTL::RR(0.25)*PI;//entropy is 12=n
      bool tmp;
      for(long a = 0;a<(1L<<n);a++)
	{
	  NTL::RR a_1 = NTL::RR(1.0);
	  NTL::RR a_2 = NTL::RR(1.0);
	  
	  long tmpa=0;
	  for(long j = 0; j<n;j++)
	    {
	      tmp =(bool)( (a&(1L<<j))>>j );
	      if(tmp)//outcome = -1
		{
		  tmpa = tmpa - (1L<<j);
		  a_1 = a_1*NTL::cos(NTL::RR(0.5)*(phi[j]+NTL::RR(0.5)*PI));
		  a_2 = (-a_2)*NTL::sin(NTL::RR(0.5)*(phi[j]+NTL::RR(0.5)*PI));
		}
	      else//outcome = +1
		{
		  tmpa = tmpa + (1L<<j);
		  a_1 = a_1*NTL::cos(NTL::RR(0.5)*(phi[j]-NTL::RR(0.5)*PI));
		  a_2 = (-a_2)*NTL::sin(NTL::RR(0.5)*(phi[j]-NTL::RR(0.5)*PI));
		}
	      
	    }
	  NTL::RR p_1 = NTL::RR(0.5)*NTL::power(a_1+a_2,2.0);
	  NTL::RR p_2 = NTL::RR(0.5)*NTL::power(a_1-a_2,2.0);
	  T.pmv.push_back( NTL::power(NTL::cos(thetasum),2.0)*p_1 + NTL::power(NTL::sin(thetasum),2.0)*p_2 );
	
	}
      delete [] theta;
      delete [] phi;
      
    }
  else if(case_nb==4)
    {
      T.nb_atoms=(1L<<10);
      NTL::RR norm_const = NTL::RR(0.0);
      
      for(long i = 3;i<T.nb_atoms+3;i++)
	{
	  if(i%(1L<<16)==0){std::cout << i << "\n";}
	  T.pmv.push_back(NTL::RR(1.0)/(NTL::RR(i)*NTL::power(NTL::log(NTL::RR(i)),1.00000001)));
	  norm_const+=T.pmv[i-3];
	  
	}
      NTL::RR tmp_ent=NTL::RR(0.0);
      
      for(long i = 3;i<T.nb_atoms+3;i++)
	{
	  
	  T.pmv[i-3]/=norm_const;
	  //std::cout << i << " --> " << T->pmv[i] << "\n";
	  
	}
      tmp_ent = tmp_ent/NTL::log(NTL::RR(2.0));
      std::cout << "Zeta-Dir related truncated prob dist entro = " << tmp_ent << "\n";
    }
  else if(case_nb==5)
    {
       //discretization of a continuous singular distribution
      
      long bit_prec = 15;//continuous distribution discretized with 2^(-bit_prec) steps size
      
      T.nb_atoms=(1L<<(bit_prec-1)); 
          
      bool tmp;
      
      NTL::RR * p = new NTL::RR [bit_prec];
      
      for(long i=0;i<bit_prec;i++)
	{
	  p[i] = NTL::RR(1.0)/(NTL::RR(1.0)+NTL::exp(-NTL::power(NTL::RR(2.0),-i)));
	  //std::cout << "\np[" << i << "] = " << p[i] ;
	}
      
      NTL::RR norm_const=NTL::RR(0.0);
      long ct=0;
      for(long i = 0 ; i < (1L<<bit_prec); i=i+2)
	{
	  T.pmv.push_back(NTL::RR(1.0));
	 
	  for(long j = 0 ; j < bit_prec ; j=j+1)//we keep even
	    {
	      tmp = (bool)((i&(1L<<j))>>j);
	      
	      if(tmp)
		{
		  T.pmv.push_back(T.pmv[ct]*(NTL::RR(1.0)-p[j]));
		  //T->remap_outcomes[i] = T->remap_outcomes[i] + 0;
		}
	      else
		{
		  T.pmv.push_back(T.pmv[ct]*p[j]);
		}
	    }
	  norm_const = norm_const + T.pmv[ct];
	  ct++;
	}
      
      for(long i = 0; i < T.nb_atoms ; i++)
	{
	  T.pmv[i] = T.pmv[i]/norm_const;
	}
      
      delete [] p;
      
    }
  else
    {
      T.nb_atoms=0;
    }
}

/**** ***** ***** ***** ***** **** ***** ***** ***** *****/

void chi2_goodness_of_fit(struct prob_dist_data T, long sample_size, long & nb_df, NTL::RR & stat)
{
  // To test is the descriptions above are correctly implemented. The sampling algorithm is agnostic
  //  to weather or not the description in an prob_dist_data is correct.
  
  // It does also timing.
  
  if(sample_size<T.nb_atoms){std::cerr << "Sample size too small for the size of the support of the distribution. --- Exit."; exit(-1);}

  long nb_classes;
  if(T.nb_atoms<=25)
    {
      nb_classes=T.nb_atoms;
    }
  else
    {
      nb_classes=25;
    }
 
  long quo = (long)std::floor(T.nb_atoms/nb_classes);//number of grouped atoms per class (maybe the last one has more, but it is the tail)
  //long rem = T.nb_atoms%nb_classes;
  if(quo==0){std::cerr << "\nDecrease number of classes. --- Exit\n\n"; exit(-1);}

  long *hist = new long [nb_classes];
  bzero(hist, nb_classes*sizeof(long));
  long index_rv;
  double time_difference=0.0;
 
  for(long i = 1; i<=sample_size;i++)
    {
      
      //std::cout << i << "\n";
      auto t_start = std::chrono::high_resolution_clock::now();
      gen_rnd_var(index_rv, T);
      auto t_end = std::chrono::high_resolution_clock::now();
      time_difference = std::chrono::duration<double, std::milli>(t_end-t_start).count();
      std::cout << i << " / " << sample_size << " : " << T.pmv[index_rv] << "\n";
      //std::cout << index_rv << " " << T.remap_outcomes[index_rv] << " " << T.pmv[index_rv] << "\n\n##########################################\n";
      for(long c = 0; c < nb_classes-1; c++)
	{
	  if( (c*quo<=index_rv) && (index_rv<(c+1)*quo) )
	    {
	      hist[c]++;
	    }
	}
      if(index_rv>=(nb_classes-1)*quo)
	{
	  hist[nb_classes-1]++;
	}
    }
  
  long sample_size_chk=0;
  for(long i = 0 ; i < nb_classes ; i++)
    {
      sample_size_chk = sample_size_chk + hist[i];
    }
  if(sample_size_chk != sample_size){std::cerr << sample_size_chk << " != " << sample_size << " problem --- exit\n"; exit(-1);} 

  NTL::RR * prob_classes = new NTL::RR [nb_classes];

  for(long i = 0; i<nb_classes;i++)
    {
      prob_classes[i]=NTL::RR(0);
    }
  for(long i = 0; i<nb_classes*quo; i++)
    {
      prob_classes[i/quo]+=T.pmv[i];
    }
  for(long i = nb_classes*quo; i<T.nb_atoms; i++)
    {
      prob_classes[nb_classes-1]+=T.pmv[i];
    }
  
  long i_left=0;
  for(long i=0;i<nb_classes;i++)
    {
      if(NTL::RR(sample_size)*prob_classes[i]>=NTL::RR(10.0))
	{
	  i_left = i;
	  break;
	}
    }
  long i_right=0;
  for(long i= nb_classes-1;i>=0;i--)
    {
      if(NTL::RR(sample_size)*prob_classes[i]>=NTL::RR(10.0))
	{
	  i_right = i;
	  break;
	}
    }

  long nb_classes_after_grouping = (i_right-i_left)+1;

  //std::cout << "\nleft index = " << i_left;
  //std::cout << "\nright index = " << i_right ;
  //std::cout << "\nnb classes grouped = " << nb_classes_after_grouping << "\n";
  
  NTL::RR * prob_classes_grouped = new NTL::RR [nb_classes_after_grouping];
  long * hist_grouped = new long [nb_classes_after_grouping];
  bzero(hist_grouped,nb_classes_after_grouping*sizeof(long));
  for(long i = 0; i<nb_classes_after_grouping;i++)
    {
      prob_classes_grouped[i]=NTL::RR(0.0);
    }

  /*
  for(long i = 0; i<nb_classes; i++)
    {
      std::cout << "\nbefore grouping --- expected freq for " << i << " class = " << NTL::RR(sample_size)*prob_classes[i];
      std::cout << "\nbefore grouping --- observed freq for " << i << " class = " << hist[i] << "\n";
      //chi2 = chi2 + NTL::power((NTL::RR(hist[i])-NTL::RR(sample_size)*prob_classes[i]),2.0)/(NTL::RR(sample_size)*prob_classes[i]);
    }
  */
  for(long i = 0; i<=i_left; i++)
    {
      hist_grouped[0]+=hist[i];
      prob_classes_grouped[0]+=prob_classes[i];
    }
 
  long itmp=1;
  for(long i = i_left+1; i<i_right; i++)
    {
      hist_grouped[itmp]+=hist[i];
      prob_classes_grouped[itmp]+=prob_classes[i];
      itmp+=1;
    }
 
  for(long i = i_right; i< nb_classes; i++)
    {
      hist_grouped[nb_classes_after_grouping-1]+=hist[i];
      prob_classes_grouped[nb_classes_after_grouping-1]+=prob_classes[i];
      
    }

  //std::cout << "\n***** ***** ***** ***** *****\n"; std::cout.flush();

  NTL::RR chi2=NTL::RR(0);
  for(long i = 0; i<nb_classes_after_grouping; i++)
     {
       /*
       std::cout << "\nafter grouping --- expected freq for " << i << " class = " << NTL::RR(sample_size)*prob_classes_grouped[i];
       std::cout << "\nafter grouping --- observed freq for " << i << " class = " << hist_grouped[i] << "\n";
       */
       chi2 = chi2 + NTL::power((NTL::RR(hist_grouped[i])-NTL::RR(sample_size)*prob_classes_grouped[i]),2.0)/(NTL::RR(sample_size)*prob_classes_grouped[i]);
     }
  nb_df = nb_classes_after_grouping-1;
  stat = chi2;

  std::cout << "\n\n***** ***** ***** ***** *****\n\n"; std::cout.flush();
  std::cout << "average time generation in milliseconds per random variable = " << time_difference/(double)sample_size << "\n";

  delete [] prob_classes_grouped;
  delete [] hist_grouped;
  delete [] prob_classes;
  delete [] hist;
}
