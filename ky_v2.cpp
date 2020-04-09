#include <iomanip>
#include <iostream>
//#include <fstream>
#include <random>
#include <ctime>
#include <chrono>
#include <cstdlib>
#include <bitset>
#include <cstring>
//#include <cassert>
//#include <unordered_map>

#include <NTL/ZZ.h>
#include <NTL/RR.h>

/**** ***** ***** ***** ***** **** ***** ***** ***** *****/

std::random_device r_dev;
std::mt19937_64 mt(r_dev());
//std::knuth_b mt(r_dev());

std::bernoulli_distribution distBer(0.5);
///std::uniform_real_distribution<double> distUnif(0, 1);
//std::uniform_int_distribution<long> distribution(0,9);

/**** ***** ***** ***** ***** **** ***** ***** ***** *****/

/*
  General rules to set the minimal accuracy.

  1) Size of the support which is the length of the probability vector, let l1 = log2(size support)
     Note here that the entropy of the probability vector <= log2(size support) so this is 
     always a safe bet whenever one does not know the exact entropy analitically.

  2) If the probability vector represents a truncation of an infinite countable discrete 
     distribution with delta as the leftover tail probilities, let l2 = -log2(-delta). Usually 
     delta is very small here.

  3) If the probability vector represents a discretization of an absolutely continuous distribution
     with 0 < epsilon < 1 as the discretization precision, let l3 = -log2(-epsilon)
     (For singular distribution, this is more a case by case analysis. Also for an absolutely continuous, 
      it is assumed that its differential entropy is bounded above so that the discretization
      is relevant.)

  4) An extra security parameter, say l4, that depends on how much resource one can use to
     sample enough in order to distinguish from the ideal (targer) distribution.

  The general rule, take min_accuray = l1 + l2 + l3 + l4;

  min_accuracy is passed as an argument to NTL::RR::SetPrecision
*/

const long min_accuracy = 1000;//ok this is a bit overkilling perhaps

const long max_number_of_coefficients=(1UL<<15);

/**** ***** ***** ***** ***** **** ***** ***** ***** *****/

/*
  To count the average number of random bits needed from
  the generator per random instance in order to compare
  with the numerically computed entropy.
*/

long ct_rnd_bit=0;

/**** ***** ***** ***** ***** **** ***** ***** ***** *****/

void roll_dice(long & outcome, long nb_faces)
{
  /*
    Outcome a uniform discrete random varible over nb_faces possible outcomes.
    It uses no more than entropy + 2 bits from the source generator according Knuth and Yao's result
    The mimimum expected number of random bits from an information point of view is log_2(nb_faces).
    Based on J. Lumbroso phd thesis idea.
    This is basically Knuth and Yao 1976 algorithm for uniform discrete distribution.
  */
  
  long x = 0;
  long y = 1;

  while(true)
    {
      y = 2*y;
      x = 2*x + distBer(mt);
      ct_rnd_bit++;
      if(y>=nb_faces)
	{
	  if(x<nb_faces)
	    {
	      outcome = x;
	      break;
	    }
	  else
	    {
	      y = y - nb_faces;
	      x = x - nb_faces;
	    }
	}
    }
}

/**** ***** ***** ***** ***** **** ***** ***** ***** *****/

struct prob_dist_data{

  long nb_atoms;//size of probability distribution support or equivalently the lenght of the probability vector
  
  std::vector<NTL::RR> pmv;//pmv.size()=nb_atoms;
  /*
    pmv = probability mass vector
    pmv[i] > 0
    sum pmv[i] = 1
    pmv[i] is accurate up to min_accuracy bits
  */

  std::vector<std::bitset<max_number_of_coefficients>> bin_rep;
  /*
    Let m = maximum of exponents
    bin_rep has this data representation:

    

                                      col 0    |  ...  |  col m-1
				      level 1  |  ...  |  level m

    row 0           (prob. mass 1) :  p_{1,1}     ...     p_{1,m}

    row 1           (prob. mass 2) :  p_{2,1}     ...     p_{2,m}

    ...

    row nb_atoms-1  (prob. mass n) :  p_{n,1}     ...     p_{n,m}
 
    ----- ----- ----- ----- -----

    bin_rep[i] = binary decomposition of pmv[i].mantissa
  */
  
  std::vector<long> nb_internal_nodes;

  std::vector<std::vector<long>> L;
  /*
    L.size = m;
    L[i].size = number of leaves at level i, type of L[i] is a hash table
    L[i][j] = k if and only if p_{k, i} = 1 for some k
    0 <= j < L[i].size
    
    nb_internal_nodes[i]= L[i].size + (number of non-leave at level i) which is
       computed with the formula. We don't want to recompute this quantity 
       for every outcome generated if we are generating a sample for instance.
       Thus it is faster to maintain this auxiliary info for generating i.i.d. sample.
  */
};

/**** ***** ***** ***** ***** **** ***** ***** ***** *****/

void entropy(NTL::RR & binary_entropy, struct prob_dist_data P)
{
  //Store the binary entropy of P in binary_entropy
  
  NTL::RR r=NTL::RR(0.0);

  for(long i = 0; i < P.nb_atoms ; i++)
    {

      r = r - P.pmv[i]*NTL::log(P.pmv[i]);
      
    }
  r = r/NTL::log(NTL::RR(2.0));
  binary_entropy = r;
}

void get_bin_rep(struct prob_dist_data & T)
{
  /*
    Get the binary representations of T.pmv and store it in T.bin_rep with T.nb_atoms rows
    and max_number_of_coefficients colums.

    If the maximum number of coefficients is not enough to hold a (min_accuracy)-bit representation/truncation,
    then the code exit and suggest the user to increase it manually. (It should seldomly occurs.)

    Get the number of internal nodes for the DDG tree reprensenting T. The representation
    of the DDG is canonical in the sens that a leaf's symbol correspond to its from left to right.
    Thus the returned outcome from sampling will be the index of a leaf. Hence
    the random outcome can be remapped.
  */
  
  if((T.nb_atoms==0)|| T.pmv.empty() ){std::cerr << "Wrong selection number --- exit\n\n"; exit(-1);}
  
  long ct_atoms=0;
  
  NTL::RR chk_sum=NTL::RR(0.0);
  
  ct_atoms=0;
  do
    {

      if(T.pmv[ct_atoms]==NTL::RR(0)){std::cerr << "zero prob --- exit --- " << ct_atoms << " " << T.pmv[ct_atoms] << "\n\n";exit(-1);}
      chk_sum += T.pmv[ct_atoms];

      std::bitset<max_number_of_coefficients> v;
      //no need to init. by default bitset field is all 0
      //A quick easy way to initialize is v = v^v;
      v = v^v;
      
      NTL::ZZ mm = T.pmv[ct_atoms].mantissa();
      long ee = T.pmv[ct_atoms].exponent();

      //check that the size of binary expansion (that is the number of coefficients)
      //  does not exceed allowed term max_number_of_coefficients. Change the latter
      //  if not enough but this should not really be the case for all practical
      //  instance in the universe. :-)
      if((unsigned long)(-ee)>=max_number_of_coefficients){std::cerr << "std::bitset length not big enough --- exit"; exit(-1);}

      long j = -ee-1;
      //std::cout << "bin rep = 0.";
      
      do
	{
	  
	  if((( mm&(NTL::ZZ(1)<<j) ) >> j) == NTL::ZZ(0))
	    {
	      v.reset(-ee-1-j);
	      //std::cout << "0";
	    }
	  else
	    {
	      v.set(-ee-1-j);
	      //std::cout << "1";
	    
	    }
	  
	  j--;
	}
      while(j>-1);
      //std::cout << "\n";

      T.bin_rep.push_back(v);
      
      ct_atoms++;
    }
  while(ct_atoms<T.nb_atoms);
  
  if( (chk_sum - NTL::RR(1.0))>= NTL::RR(T.nb_atoms)*NTL::power(NTL::RR(2.0),-min_accuracy) ) {std::cerr << "check sum prob vector error --- exit\n\n"; exit(-1);}
  if(T.bin_rep.size() != (unsigned long)(ct_atoms) ) {std::cerr << "big shit --- exit\n\n"; exit(-1);}

  /****************************************/
  
  for(unsigned long l=0;l<max_number_of_coefficients;l++)
    {
      std::vector<long> v;
      
      for(long i = 0; i<T.nb_atoms;i++)
	{
	  if(T.bin_rep[i][l]==true)
	    {
	      v.push_back(i);
	    }
	}
      T.L.push_back(v);
      
    }

  /****************************************/
  
  T.nb_internal_nodes.push_back(2);//(there are always 2 nodes for the 1st level indexed by (1 - 1) )
  for(long l=1;l<min_accuracy;l++)//only the first min_accuracy bits out of the max_number_of_coefficents are accurate, do not loop until max_number_of_coefficients!
    {
      //T->nb_nodes[l-1]=(1L<<l);//maximum possibles of nodes at level l indexed by (l-1) + 1, 
      long tmp = 0;
      for(long m=0;m<l;m++)//account of previous level
	{
	  tmp = tmp + (1L<<(l-m))*(T.L[m].size());
	}
      T.nb_internal_nodes.push_back((1L<<(l+1))-tmp);
    }
}

/**** ***** ***** ***** ***** **** ***** ***** ***** *****/

void gen_rnd_var(long & index_rv, struct prob_dist_data T)
{
   /*
    index_rv is a random outcome from distribution described by T.
    0 <= index_rv < size of support of T = T.nb_atoms

    A user can transform the range of index_rv if required. For instance by an affine transformation.
    
    It uses no more than the entropy(T) + 2 bits on average from the random source.
    A simple adaptation of the roll_dice above using lists.
    This is basically Knuth and Yao 1976 algorithm for discrete distributions.
   */
  
  long x=0;
  long y=1;
  long lev=0;
  
  while(true)
    {
      x = 2*x + distBer(mt);
      ct_rnd_bit++;
      
      y = 2*y;
      if( y >= T.nb_internal_nodes[lev])//number of internal nodes = exit nodes (leaves) + non-exit nodes (continue to next level)
	{
	  
	  if(x<(long)T.L[lev].size())//number of leaves
	    {
	      
	      index_rv = T.L[lev][x];
	      return;
	    
	    }
	  else
	    {
	      y = y - T.L[lev].size();
	      x = x - T.L[lev].size();
	    }
	}
      lev++;
    }
}

/**** ***** ***** ***** ***** **** ***** ***** ***** *****/

/*
  Descriptions of a few probability distributions together with a chi square test.

  Generic sampling algorithms like above are generally agnostic to the description
    of the sampled distribution. Therefore the chi sqaure test can be use to validate
    the implementation (assuming the chi sqaure test is correctly implemented.)

    Add any of your probability distributions in prob_dist_list.h
*/

#include "prob_dist_list.h" 
#include "ky_v2_IO.h"

int main(void)
{
 
  std::cout.precision(15);
  std::cout.setf( std::ios::fixed, std::ios::floatfield );

  NTL::RR::SetPrecision(min_accuracy);
  NTL::RR::SetOutputPrecision(15);

  struct prob_dist_data P;

  //See prob_dis_list.h for descriptions
  //which_prob_dist(1,P);

  query_user(P);
  
  get_bin_rep(P);

  /*
    Example of a call:
  */
  if(false)
    {
      long rv;
      gen_rnd_var(rv,P);
      
      std::cout << "\nP{outcome = "<<rv<<"} = " << P.pmv[rv] << "\n";
    }
  
  long nb_df;//degrees of freedom resulting from chi square test
  NTL::RR stat_test_value;

  /*
    Make it big enough to make sure the implementation is checked properly. 
    The chi sqaure groups classes if the observed frequence < 10, and hence
    resulting in a decrease of the number of degrees of freedom.
    The number of classes is set to 25, that it splits number of possible outcomes 
    into 25 classes.
  */
  long sample_size;
  if(P.nb_atoms < 25*10)
    {
      sample_size = 250;
    }
  else
    {
      sample_size=P.nb_atoms*10;
    }
  
  chi2_goodness_of_fit(P, sample_size, nb_df, stat_test_value);

  std::cout << "\n\nnumber of atoms = " << P.nb_atoms;
  
  NTL::RR bin_ent;
  entropy(bin_ent,P);
  std::cout << "\nentropy = " << bin_ent;
  
  std::cout << "\n\nsample size = " << sample_size;
  
  std::cout.flush();
  
  std::cout << "\naverage number of random bits used per random variable = " << (double)ct_rnd_bit/(double)sample_size;
  std::cout << " (compare with entropy)";
  std::cout << "\ngoodness of fit chi2 stat (" << nb_df << " d.f.) = " << stat_test_value << "\n\n" ;
    
  return 0;
}
