void query_user(struct prob_dist_data & T)
{
  
  std::string filename;
  std::string choice;

  std::cout <<
	 "Choose among:\n" <<
	 "  (1) arbitrary mass vector\n" <<
	 "  (2) discrete gaussian\n" <<
	 "  (3) binomial\n" <<
	 "  (4) Zeta-Dirichlet\n" <<
	 "  (5) GHZ quantum mechanical\n" <<
         "  (6) Discretization of a singular continuous distribution\n" <<
	 "  anything else --- exit\n";
  
  std::cout << "\nEnter choice: ";
  std::getline(std::cin , choice);

  bool chk_valid_choice=false;
  for(long i = 1 ; i <= 6 ; i++)
    {
      chk_valid_choice = chk_valid_choice || (std::stol(choice)!=i) ;
    }
  if(!chk_valid_choice)
    {
      std::cerr << "wrong choice --- exit\n";
    }
  
  if(std::stol(choice) == 1)
    {
      std::string tmpstr;
      std::cout << "\n\nSize of the support: ";
      
      getline(std::cin, tmpstr);
      if(std::stod(tmpstr)<=1)
	{
	  std::cerr << "error --- support size <=1 \n";
	  exit(-1);
	}
      
      T.nb_atoms = std::stod(tmpstr);
      
      for(long i = 0 ; i < T.nb_atoms; i++)
	{
	  std::cout << "\n Prob{Outcome "<<i<<"} = ";
	  getline(std::cin, tmpstr);
	  T.pmv.push_back(NTL::RR(std::stod(tmpstr)));
	}
      std::cout << "\n";
      std::cout.flush();
    }
  else if(std::stol(choice) == 2)
    {
      std::string tmpstr;
      std::cout << "\n\nAverage = ";
      getline(std::cin,tmpstr);
      double alp = std::stod(tmpstr);

      std::cout << "\nStandard deviation = ";
      getline(std::cin,tmpstr);
      double bet = std::stod(tmpstr);
      if(bet<=0)
	{
	  std::cerr << "\nstardard deviation <= 0 --- exit\n";
	  exit(-1);
	}

      std::cout << "\nHow far away, in a number of standard deviation, from the mean are the truncation points? ";
      getline(std::cin,tmpstr);
      double dis = std::stod(tmpstr);
      if(dis<=0)
	{
	  std::cerr << "\ndistance <= 0 --- exit\n";
	  exit(-1);
	}
      
      NTL::RR alpha = NTL::exp(NTL::RR(alp));
      NTL::RR beta  = NTL::RR(bet);
      
      NTL::ZZ zleft = NTL::conv<NTL::ZZ>(NTL::floor(alpha - NTL::RR(dis)*beta));
      NTL::ZZ zright= NTL::conv<NTL::ZZ>(NTL::ceil(alpha + NTL::RR(dis)*beta));
      
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
  else if(std::stol(choice)==3)
    {
      std::string tmpstr;
      std::cout << "\n\nNumber of trials = ";
      getline(std::cin,tmpstr);
      long NB_TRIALS = std::stol(tmpstr);
      if(NB_TRIALS<=1)
	{
	  std::cerr << "\nnumber of trials <= 1 --- exit\n\n";
	  exit(-1);
	}

      std::cout << "\nSuccess probability = ";
      getline(std::cin, tmpstr);
      double sp = std::stod(tmpstr);
      if(sp < 0 || sp > 1)
	{
	  std::cerr << "\nsuccess probability <= 0 or >= 1 --- exit\n\n";
	  exit(-1);
	}
      
   
      NTL::RR SUCCESS_PROB = NTL::RR(sp);

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
  else if(std::stol(choice)==4)
    {

      std::cout << "\n\nConcentration parameter = ";
      std::string tmpstr;
      getline(std::cin,tmpstr);
      double u = std::stol(tmpstr);
      if(u<=0.0)
	{
	  std::cout << "\n\nconcentration parameter <= 0 --- exit\n";
	  exit(-1);
	}
      std::cout << "\n\nTruncate where ? ";
      getline(std::cin,tmpstr);
      long nbat = std::stol(tmpstr);
      if(nbat<=1)
	{
	  std::cout << "\n\ncut off point <= 1 --- exit\n";
	  exit(-1);
	}
      
      T.nb_atoms=nbat;
      NTL::RR norm_const = NTL::RR(0.0);
      
      for(long i = 3;i<T.nb_atoms+3;i++)
	{
	  if(i%(1L<<16)==0){std::cout << i << "\n";}
	  T.pmv.push_back(NTL::RR(1.0)/(NTL::RR(i)*NTL::power(NTL::log(NTL::RR(i)),u)));
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
  else if(std::stod(choice)==5)
    {
      std::cout << "\n\nInput n where Hilbert space dimension = 2^n : ";
      std::string tmpstr;
      getline(std::cin,tmpstr);
      
      long n = std::stol(tmpstr);
      if(n<1)
	{
	  std::cerr << "\n\nn < 1 --- exit\n\n";
	  exit(-1);
	}
      
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
  else if(std::stol(choice)==6)
    {
      //discretization of a continuous singular distribution
      std::cout << "\n\nInput p>1 for 2^(-p) = step size of discretization: ";
      std::string tmpstr;
      getline(std::cin,tmpstr);
      long bit_prec = std::stol(tmpstr);//continuous distribution discretized with 2^(-bit_prec) steps size
      if(bit_prec<=1)
	{
	  std::cerr << "\np <= 1 --- exit\n\n";
	  exit(-1);
	}
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
      std::cerr << "\n\nsmelly onara --- sayonara\n\n";
      exit(-1);
    }
}




