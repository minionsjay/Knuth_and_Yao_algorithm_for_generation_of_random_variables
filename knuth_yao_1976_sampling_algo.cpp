/***** BEGINNING OF DISCLAIMER *****
       
       This code is provided without warranty of any kind, either express or implied. Use at your own risk.
       
       The use of this is done at your own discretion and risk and with agreement that you will be solely responsible for any damage to your computer system or loss of data that results from the use of this code. You are solely responsible for adequate protection and backup of the data and equipment used in connection with this code, and we will not be liable for any damages that you may suffer in connection with using, modifying or distributing any of this code. No advice or information, whether oral or written, obtained by you from us or from this website shall create any warranty for this code.
       
       We make makes no warranty that:
       
       1) the code will meet your requirements
       2) the code will be uninterrupted, timely, secure or error-free
       3) the results that may be obtained from the use of the software will be effective, accurate or reliable
       4) the quality of the code will meet your expectations
       5) any errors in the code obtained from us will be corrected. 
       
       The code and its documentation made available on this website:
       
       1) could include technical or other mistakes, inaccuracies or typographical errors. We may make changes to the code or documentation made available on its web site at any time without prior-notice.
       2) may be out of date, and we make no commitment to update such materials. 
       
       We assume no responsibility for errors or omissions in the code or documentation available from its web site.
       
       In no event shall we be liable to you or any third parties for any special, punitive, incidental, indirect or consequential damages of any kind, or any damages whatsoever, including, without limitation, those resulting from loss of use, data or profits, and on any theory of liability, arising out of or in connection with the use of this code. 
       
       ***** END OF DISCLAIMER *****/

/*
  
  AUTHOR: Claude Gravel
  
  VERSION: 1 (March 2020)
  
  TECHNICAL REFERENCE: ky_algo.pdf

  COMPILE AND LINK COMMAND EXAMPLE: "g++ -Wall -Wpedantic -g -O5 knuth_yao_1976_sampling_algo.cpp -o knuth_yao_1976_sampling_algo -std=c++17"

  RUN COMMAND EXAMPLE: "./knuth_yao_1976_sampling_algo" 
  
*/

#include <iomanip>
#include <iostream>
#include <fstream>
#include <random>
#include <time.h>
#include <sys/resource.h>

#include <cstdlib>
#include <cstring>
#include <cassert>

/**** ***** ***** ***** *****/

std::random_device r_dev;//source of randomness
std::mt19937_64 mt(r_dev());//type of pseudo-random bit generator
//std::knuth_b mt(r_dev());

std::bernoulli_distribution distBer(0.5);//generate unbiased Bernoulli i.i.d. random variales

/**** ***** ***** ***** *****/

struct prob_dist_basic_dat
{
  long nb_atoms;//size of the support of the probability distribution that is the lenght of the probability vector
  long min_accuracy;//minimal accuracy guaranteed by the user for the probabilities
  long max_size_storage_bin_rep;//must be power of 2 for storage of rows (binary representation) of brt below
  //max_size_storage_bin_rep should be sufficiently large to hold min_accuracy bits
  
  bool ** brt;//table for binary representation, brt[i][j] = j-th coefficient of binary representation of i-th probability value
  long * nb_nodes;//nb_nodes[j] = number of possible nodes at j-th level, j-th level is indexed by (j-1) 
  long * p2_val_nb_nodes;//ceiling log base two of number of nodes
  long * nb_exit_nodes;//p2_val_nb_nodes[j] + nb_exit_nodes[j] = nb_nodes[j]
};

/**** ***** ***** ***** *****/

long ct_rnd_bit=0;//just to count how many time distBer is called and then averaged over the sample where sample size is specified by the user together with the file containing the proper binary representations. The average should between H and H + 2, where H is the entropy of the distribution.

/**** ***** ***** ***** *****/

double get_ms_res_time(void)//Added here for timing purpose, the higher the entropy is, the longer it takes to sample for a given size of the latter.
{
  
  double current_time;
  
  struct rusage ruse;
  getrusage(RUSAGE_THREAD, &ruse);
  current_time = ((double) ruse.ru_utime.tv_sec * 1000.0) + (((double) ruse.ru_utime.tv_usec) / 1000.0); //resultat en millisecondes
  
  return current_time;
}

/**** ***** ***** ***** *****/

void roll_dice(long * outcome, long nb_faces)//Jeremy Lumbroso adaptation of Knuth-Yao algorithm for the uniform distribution. This procedure is added here for reference purpose, but is not needed anywhere else. 
{
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
	      *outcome = x;
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

/**** ***** ***** ***** *****/

void read_bin_rep_from_file(std::string fname, struct prob_dist_basic_dat * T)
{
  /*
    Read data representing a probability distribution from a binary file.

    The format is simple and as follow:
    
    Bytes 0  to  7 inclusively = long integer encoding the number of atoms of the distribution that is the length
                               = T->nb_atoms;

    Bytes 8  to 15 inclusively = minimal accuracy of the probability values
                               = T->min_accuracy

    Bytes 16 to 23 inclusively = maximal length in bits of binary representation and must be a power of two > minimal accuracy
                               = T->max_size_storage_bin_rep
    
    Bytes 24 + i*(T->max_size_storage_bin_rep) to 24 + (i+1)*(T->max_size_storage_bin_rep) - 1 inclusively = binary representation of the i-th probability value for 0 <= i < T->nb_atoms

   */
  
  if(fname.empty()){std::cerr << "problem --- exit\n"; exit(-1);}
 
  std::ifstream file_res(fname, std::ifstream::in | std::ifstream::binary);
  assert((file_res));
  
  file_res.read((char *)(&T->nb_atoms),sizeof(long));
  file_res.read((char *)(&T->min_accuracy),sizeof(long));
  file_res.read((char *)(&T->max_size_storage_bin_rep),sizeof(long));
  T->brt = new bool * [T->nb_atoms];
  
  for(long i=0;i<T->nb_atoms;i++)
    {
      T->brt[i] = new bool [T->max_size_storage_bin_rep];
      file_res.read((char *)((T->brt[i])),(T->max_size_storage_bin_rep/8));
    }
  file_res.close();
  
  /****************************************/
  
  T->nb_exit_nodes = new long [T->min_accuracy];
  bzero(T->nb_exit_nodes,(T->min_accuracy)*sizeof(long));
  for(long l=0;l<T->min_accuracy;l++)
    {
      for(long i = 0; i<T->nb_atoms;i++)
	{
	  T->nb_exit_nodes[l]=T->nb_exit_nodes[l]+T->brt[i][l];
	}
    }
  
  /****************************************/
  
  T->nb_nodes = new long [T->min_accuracy];
  bzero(T->nb_nodes,(T->min_accuracy)*sizeof(long));
  
  T->nb_nodes[0]=2;//(there are always 2 nodes for the 1st level indexed by (1 - 1) )
  for(long l=1;l<T->min_accuracy;l++)
    {
      //T->nb_nodes[l-1]=(1L<<l);//maximum possibles of nodes at level l indexed by (l-1) + 1, 
      long tmp = 0;
      for(long m=0;m<l;m++)//account of previous level
	{
	  //T->nb_nodes[l-1]=T->nb_nodes[l-1] - ( (1L<<(l-m))*(T->nb_exit_nodes[l-m-1]) );
	  tmp = tmp + (1L<<(l-m))*(T->nb_exit_nodes[m]);
	}
      T->nb_nodes[l]=(1L<<(l+1))-tmp;
    }
  
  T->p2_val_nb_nodes = new long [T->min_accuracy];
  bzero(T->p2_val_nb_nodes,(T->min_accuracy)*sizeof(long));
  
  for(long l=0;l<T->min_accuracy;l++)
    {
      T->p2_val_nb_nodes[l]=(long)std::ceil(std::log2(T->nb_nodes[l]));
    }
  
}

/**** ***** ***** ***** *****/

void gen_rnd_var_KY(long * index_rv, struct prob_dist_basic_dat T)
{
  /*
    This is the Knuth and Yao 1976 algorithm in procedural form.
    
    index_rv is a pointer to the generated random variable.
    
    T is the description of the probability distribution through the data structure prob_dist_basic_dat.
  */
  
  long x=0;
  long y=1;
  long lev=0;
  long tmp = 0;
  bool chk_out=false;
  while(!chk_out)
    {
      x = 2*x + distBer(mt);
      ct_rnd_bit++;
      
      y = 2*y;
      if( y >= T.nb_nodes[lev])
	{
	  
	  if(x<T.nb_exit_nodes[lev])
	    {
	      tmp = 0;
	      for(long k = 0; k < T.nb_atoms ; k++)
		{
		  tmp = tmp + T.brt[k][lev];
		  if(x==tmp-1)
		    {
		      *index_rv=k;
		      chk_out=true;
		      break;
		    }
		}
	    }
	  else
	    {
	      y = y - T.nb_exit_nodes[lev];
	      x = x - T.nb_exit_nodes[lev];
	    }
	}
      lev++;
    }
}

/**** ***** ***** ***** *****/

int main(void)
{
  std::cout.precision(15);
  std::cout.setf( std::ios::fixed, std::ios::floatfield );

  /****************************************/
  
  struct prob_dist_basic_dat Q;

  std::string filename;
  std::cout << "Enter file name: ";
  getline(std::cin , filename);
  if(!filename.empty())
    {
      read_bin_rep_from_file(filename, &Q);      
    }
  else
    {
      std::cerr << "empty file name --- exit\n\n";
      exit(-1);
    }
  
  long rv;
  long sample_size;
  std::cout << "\nEnter sample size: ";
  std::cin >> sample_size;
  if(sample_size<1)
    {
      std::cerr << "sample size < 1 --- exit\n\n";
      exit(-1);
    }
  double t0=0.0;
  double t1=0.0;
  for(long i=0;i<sample_size;i++)
    {
      t0 += get_ms_res_time();
      gen_rnd_var_KY(&rv, Q);
      t1 += get_ms_res_time();
      //std::cout << i << " " << rv << "\n";
    }

  std::cout << "\n\nLength of the probability vector = " << Q.nb_atoms;
  std::cout << "\nMinimal binary accuracy per probabity value = " << Q.min_accuracy;

  std::cout << "\n\nEstimation of entropy = " << (double)ct_rnd_bit/(double)sample_size ;
  std::cout << "\n(Note that for not too large sample, the estimation should be between H and H + 2\nwith very high confidence where H is the entropy.)" ;
  std::cout << "\n\nAverage time generation per random variable = " << (t1-t0)/(double)sample_size << "\n\n";

  for(long i = 0;i<Q.nb_atoms;i++)
    {
      delete [] Q.brt[i];
    }
  delete [] Q.brt;

  delete [] Q.nb_nodes;
  delete [] Q.p2_val_nb_nodes;
  delete [] Q.nb_exit_nodes;
  
  return 0;
}
