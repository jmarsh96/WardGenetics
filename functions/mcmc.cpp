#include <Rcpp.h>
#include <iostream>
#include <fstream>
using namespace Rcpp;

// [[Rcpp::plugins("cpp11")]]

IntegerVector WhichVec(int x, IntegerVector vec);
int ReturnSourceSequence(int sequence_loc, List data);
IntegerVector ReturnDuplicateSequences(int sequence_loc, List data);
IntegerVector ReturnGenSourceVector(List data);
List UpdateImputedNodesMoveAll(List data_cur, List data_can, NumericVector parameters);
List ReturnNodesToImpute(List data);
IntegerVector RemoveElement(int loc, IntegerVector vec);

List InitialiseData(IntegerVector t_a, IntegerVector t_c, IntegerVector t_d, IntegerVector source, IntegerVector hcw_ind, IntegerMatrix screening_matrix)
{
  IntegerVector zero_based_source = clone(source);
  for(int i = 0; i<source.length(); i++)
  {
    if(source[i] > 0)
    {
      zero_based_source[i]--;
    }
  }
  
  int num_patients = t_c.length();
  IntegerVector hcw_loc = WhichVec(1, hcw_ind);
  if(hcw_loc.length() > 0) {
    num_patients = hcw_loc[0];
  } 
  
  List data = List::create(Named("t_a") = t_a , 
                           Named("t_c") = t_c,
                           Named("t_d") = t_d,
                           Named("source") = zero_based_source,
                           Named("screening_matrix") = screening_matrix,
                           Named("num_patients") = num_patients,
                           Named("hcw_ind") = hcw_ind);
  return data;
}


// [[Rcpp::export]]
List InitialiseData_GEN(int sequence_length, IntegerVector t_a, IntegerVector t_c, IntegerVector t_d, IntegerVector source, IntegerVector hcw_ind, 
                        IntegerMatrix screening_matrix, IntegerVector genetic_ids, IntegerVector sample_times, IntegerVector variant_numbers, 
                        IntegerMatrix genetic_matrix)
{
  Rcout << "good" << std::endl;
  IntegerVector zero_based_source = clone(source);
  for(int i = 0; i<source.length(); i++)
  {
    if(source[i] > 0)
    {
      zero_based_source[i]--;
    }
  }
  IntegerVector zero_based_genetic_id = clone(genetic_ids);
  for(int i = 0; i<genetic_ids.length(); i++)
  {
    zero_based_genetic_id[i]--;
  }
  
  int num_patients = t_c.length();
  IntegerVector hcw_loc = WhichVec(1, hcw_ind);
  if(hcw_loc.length() > 0) {
    num_patients = hcw_loc[0];
  } 
  
  
  IntegerVector imputed_nodes(genetic_ids.length());
  List data = List::create(Named("t_a") = t_a , 
                           Named("t_c") = t_c,
                           Named("t_d") = t_d,
                           Named("source") = zero_based_source,
                           Named("screening_matrix") = screening_matrix,
                           Named("genetic_ids") = zero_based_genetic_id,
                           Named("sample_times") = sample_times,
                           Named("variant_numbers") = variant_numbers,
                           Named("genetic_matrix") = genetic_matrix,
                           Named("sequence_length") = sequence_length,
                           Named("imputed_nodes") = imputed_nodes,
                           Named("num_patients") = num_patients,
                           Named("hcw_ind") = hcw_ind);
  return data;
}


IntegerVector WhichVec(int x, IntegerVector vec)
{
  IntegerVector locs;
  for(int i = 0; i<vec.length(); i++)
  {
    if(vec[i]==x)
    {
      locs.push_back(i);
    }
  }
  return locs;
}

IntegerVector WhichNotEqualVec(int x, IntegerVector vec)
{
  IntegerVector locs;
  for(int i = 0; i<vec.length(); i++)
  {
    if(vec[i]!=x)
    {
      locs.push_back(i);
    }
  }
  return locs;
}

IntegerVector RemoveElement(int loc, IntegerVector vec)
{
  IntegerVector new_vec;
  for(int i = 0; i<vec.length(); i++)
  {
    if(i!=loc) new_vec.push_back(vec[i]);
  }
  return new_vec;
}

int SampleVector(IntegerVector x)
{
  double U = R::runif(0.0,1.0);
  int length = x.length();
  return x[floor(U*length)];
}

// [[Rcpp::export]]
IntegerVector SampleMultipleVector(IntegerVector x, int size)
{
  if(size > x.length()) stop("Too few elements to sample");
  
  IntegerVector out;
  for(int i = 0; i<size; i++)
  {
    double U = R::runif(0.0,1.0);
    int length = x.length();
    out.push_back(x[floor(U*length)]);
    x = RemoveElement(floor(U*length), x);
  }
  return out;
}



int col_pop(List data, int day)
{
  int cp = 0;
  int num_patients = data["num_patients"];
  IntegerVector t_c = data["t_c"];
  IntegerVector t_d = data["t_d"];
  IntegerVector source = data["source"];
  for(int i = 0; i<num_patients; i++)
  {
    //Rcout << "Day = " << day << ", i = " << i << ", t_c = " << t_c[i] << ", t_d"
    if(t_c[i] <= day && t_c[i] != -1 && t_d[i] >= day)
    {
      // The individual has a colonisation time on or before the day and are still present in the ward
      if(t_c[i]==day && source[i] == -1)
      {
        // The individual was imported on the day
        cp++;
      }
      else if(t_c[i] < day)
      {
        // The individual was colonised prior to the day
        cp++;
      }
    }
  }
  return cp;
}

void PrintNumVectorToFile(std::ofstream &myfile, NumericVector vec)
{
  for(int i = 0; i<vec.length(); i++)
  {
    myfile << vec[i] << " ";
  }
  myfile << std::endl;
}

void PrintIntVectorToFile(std::ofstream &myfile, IntegerVector vec)
{
  for(int i = 0; i<vec.length(); i++)
  {
    myfile << vec[i] << " ";
  }
  myfile << std::endl;
}

// Runs a quick check that there are no negative values in the distance matrix
bool DoesMatrixContainNegativeEntries(IntegerMatrix mat)
{
  for(int i = 0; i<mat.nrow(); i++)
  {
    for(int j = 0; j<mat.ncol(); j++)
    {
      if(mat(i,j)<0) return true;
    }
  }
  return false;
}

// [[Rcpp::export]]
IntegerVector ReturnPathToRoot_Rcpp(IntegerVector gen_source, int node)
{
  int max_length = gen_source.length();
  int counter = 1;
  IntegerVector path;
  int parent = gen_source[node];
  path.push_back(node);
  path.push_back(parent);
  bool path_found = false;
  if(parent == -1)
  {
    return path;
  }
  while(!path_found)
  {
    node = parent;
    parent = gen_source[node];
    path.push_back(parent);
    if(parent == -1)
    {
      path_found = true;
    }
    counter++;
    if(counter > max_length)
    {
      Rcout << "Problem node = " << node << std::endl;
      Rcout << "Gen source = " << gen_source << std::endl;
      stop("Problem with the genetic source vector, there is a loop!");
    }
  }
  return path;
}

int CalculateDistanceBetweenNodes_Rcpp(int node1, int node2, IntegerVector gen_source, IntegerMatrix genetic_matrix)
{
  if(node1 == node2) return 0;
  bool print_debug = false;
  IntegerVector path1 = ReturnPathToRoot_Rcpp(gen_source, node1);
  IntegerVector path2 = ReturnPathToRoot_Rcpp(gen_source, node2);
  int distance = 0;
  // Assign path1 to be the longest path
  if(path1.length() < path2.length())
  {
    IntegerVector temp_path = path1;
    path1 = path2;
    path2 = temp_path;
    
    int temp_node = node1;
    node1 = node2;
    node2 = temp_node;
  }
  if(print_debug)
  {
    Rcout << "path 1 = " << path1 << std::endl;
    Rcout << "path 2 = " << path2 << std::endl;
  }
  
  // path 1 is longer or equal to path 2
  // Check if node 2 is contained in path 1
  IntegerVector node2_path1_loc = WhichVec(node2, path1);
  if(node2_path1_loc.length() != 0)
  {
    // Node two is contained within path one, sum distances
    int parent_node = gen_source[node1];
    distance += genetic_matrix(node1, parent_node);
    while(parent_node != node2)
    {
      int current_node = parent_node;
      parent_node = gen_source[current_node];
      distance += genetic_matrix(current_node, parent_node);
    }
    return distance;
  }
  
  // they are not contained in the path, therefore try to find the node of divergence
  int path1_loc, path2_loc;
  bool common_node_found = false;
  for(int i = 0; i<path2.length(); i++)
  {
    int current_node = path2[i];
    IntegerVector cur_node_loc_in_path1 = WhichVec(current_node, path1);
    if(cur_node_loc_in_path1.length() != 0)
    {
      path1_loc = cur_node_loc_in_path1[0];
      path2_loc = i;
      common_node_found = true;
      if(print_debug) Rcout << "path1 loc = " << path1_loc << ", path2 loc = " << path2_loc << std::endl;
      break;
    }
  }
  
  
  if(common_node_found)
  {
    int common_node = path2[path2_loc];
    if(print_debug) Rcout << "Common node = " << common_node << std::endl;
    if(common_node == -1)
    {
      path1_loc--;
      path2_loc--;
      distance += genetic_matrix(path1[path1_loc], path2[path2_loc]);
      if(print_debug) Rcout << "Adding nodes " << path1[path1_loc] << " and " << path2[path2_loc] << std::endl;
    }
    
    for(int i = 0; i<path1_loc; i++)
    {
      distance += genetic_matrix(path1[i], path1[i+1]);
      if(print_debug) Rcout << "i = " << i << ", " << path1[i] << " -> " << path1[i+1] << ", dist = " << genetic_matrix(path1[i], path1[i+1]) << std::endl;
    }
    
    for(int i = 0; i<path2_loc; i++)
    {
      distance += genetic_matrix(path2[i], path2[i+1]);
      if(print_debug) Rcout << "i = " << i << ", " << path1[i] << " -> " << path1[i+1] << ", dist = " << genetic_matrix(path1[i], path1[i+1]) << std::endl;
      
    }
    return distance;
  }
  else
  {
    for(int i = 0; i<path1.length()-1; i++)
    {
      distance += genetic_matrix(path1[i], path1[i+1]);
    }
    
    for(int i = 0; i<path2.length()-1; i++)
    {
      distance += genetic_matrix(path2[i], path2[i+1]);
    }
    return distance;
  }
}

void PrintColonisationTimeSumToFile(std::ofstream &myfile, List data)
{
  IntegerVector t_c = data["t_c"];
  int colonisation_time_sum = 0;
  for(int i = 0; i<t_c.length(); i++)
  {
    if(t_c[i] != -1)
    {
      colonisation_time_sum += t_c[i];
    }
  }
  myfile << colonisation_time_sum << std::endl;
}

// [[Rcpp::export]]
int hcw_col_pop(List data, int day)
{
  int cp = 0;
  int num_patients = data["num_patients"];
  IntegerVector t_c = data["t_c"];
  IntegerVector t_d = data["t_d"];
  IntegerVector source = data["source"];
  
  // Quick check to exit early if there are no hcw present
  //Rcout << std::endl << "num patients = " << num_patients << ", tc length = " << t_c.length() << std::endl;
  if(num_patients == t_c.length()) return 0;
  
  for(int i = num_patients; i<t_c.length(); i++)
  {
    //Rcout << "Day = " << day << ", i = " << i << ", t_c = " << t_c[i] << ", t_d" << t_d[i] << std::endl;
    if(t_c[i] < day && t_c[i] != -1 && t_d[i] >= day)
    {
      cp++;
    }
  }
  return cp;
}

int total_col_pop(List data, int day)
{
  return hcw_col_pop(data, day) + col_pop(data, day);
}

int CalculateTruePositives(List data)
{
  IntegerMatrix screening_matrix = data["screening_matrix"];
  int num_patients = data["num_patients"];
  int TP = 0;
  for(int i = 0; i<num_patients; i++)
  {
    IntegerVector y = WhichVec(1,screening_matrix(i,_));
    TP += y.length();
  }
  return TP;
}

// [[Rcpp::export]]
int CalculateFalseNegatives(List data)
{
  int num_patients = data["num_patients"];
  IntegerVector t_c = data["t_c"];
  IntegerVector t_d = data["t_d"];
  IntegerMatrix screening_matrix = data["screening_matrix"];
  int FN = 0;
  for(int i = 0; i<num_patients; i++)
  {
    if(t_c[i] != -1)
    {
      for(int t = t_c[i]; t<=t_d[i]; t++)
      {
        if(screening_matrix(i,t)==0)
        {
          FN++;
        }
      }
    }
  }
  return FN;
}

double CalculateScreeningLikelihood(List data, NumericVector parameters)
{
  double z = parameters[0];
  double loglik = CalculateTruePositives(data) * log(z) + CalculateFalseNegatives(data) * log(1-z);
  //Rcout << "TP = " << CalculateTruePositives(data) << ", FN = " << CalculateFalseNegatives(data) << std::endl;
  return loglik;
}

double CalculateImportationLikelihood(List data, NumericVector parameters)
{
  double p = parameters[1];
  IntegerVector source = data["source"];
  IntegerVector import_idx = WhichVec(-1,source);
  int population_size = data["num_patients"];
  int importation_sum = import_idx.length();
  double loglik = importation_sum * log(p) + (population_size-importation_sum) * log(1-p);
  //Rcout << "Import sum = " << importation_sum << ", p = " << p << ", pop size = " << population_size << std::endl;
  return loglik;
}

double CalculateMutationProbability(double rate, double t)
{
  double prob = 0.75*(1-exp(-4*rate*t));
  return prob;
}



double CalculateGeneticLikelihood(List data, NumericVector parameters)
{
  double lambda = parameters[4];
  double importation_distance = parameters[5];
  //double variant_distance = parameters[5];
  //Rcout << "Try calculate gen source " << std::endl;
  IntegerVector gen_source = ReturnGenSourceVector(data);
  //Rcout << "In Gen ll - Gen source = " << gen_source << std::endl;
  IntegerVector sample_times = data["sample_times"];
  IntegerMatrix genetic_matrix = data["genetic_matrix"];
  int sequence_length = data["sequence_length"];
  if(gen_source.length() != genetic_matrix.nrow()) stop("Inconsistent length of genetic ids/matrix");
  double loglik = 0.0;
  int zero_time_counter = 0;
  for(int i = 0; i < gen_source.length(); i++)
  {
    int cur_source = gen_source[i];
    if(cur_source >= 0)
    {
      // Check they have a source and are not an importation
      int dist = genetic_matrix(cur_source,i);
      int time_diff = sample_times[i] - sample_times[cur_source];
      //Rcout << "i = " << i << ", source = " << cur_source << ", dist = " << dist << ", time diff = " << time_diff << std::endl;
      if(time_diff > 0)
      {
        double mutation_probability = CalculateMutationProbability(lambda, time_diff);
        loglik += (sequence_length-dist)*log(1-mutation_probability) + dist*log(mutation_probability/3);
      }
      else
      {
        zero_time_counter++;
        //double mutation_probability = CalculateMutationProbability(lambda, 0.01);
        //loglik += (sequence_length-dist)*log(1-mutation_probability) + dist*log(mutation_probability/3);
        if(dist != 0)
        {
          //Rcout << "Non zero distance from target " << i << " and genetic source " << cur_source
          //      << ", with time difference of " << time_diff << std::endl;
          loglik -= 999999999999;
        }
      }
    }
    else if(cur_source == -1)
    {
      // the sequence does not have a previous source, therefore look to the vector of imputed distances
      //int dist = master_distances[i];
      //loglik += R::dbinom(dist, sequence_length, importation_p, 1);
    }
  }
  IntegerVector variant_numbers = data["variant_numbers"];
  IntegerVector import_idx = WhichVec(-1, gen_source);
  //IntegerVector import_variant_numbers = as<IntegerVector>(variant_numbers[import_idx]);
  //Rcout << "Import idx = " << import_idx << std::endl;
  if(import_idx.length() > 1)
  {
    for(int i = 0; i<(import_idx.length()-1); i++)
    {
      for(int j = (i+1); j<import_idx.length(); j++)
      {
        int ii = import_idx[i];
        int jj = import_idx[j];
        int dist = genetic_matrix(ii,jj);
        //Rcout << "i = " << i << ", j = " << j << ", dist = " << dist << std::endl;
        double rate;
        if(variant_numbers[ii]==variant_numbers[jj])
        {
          rate = importation_distance;
          if(rate > 0)
          {
            loglik += R::dpois(dist, rate, 1);
          }
          
        }
        else
        {
          //rate = variant_distance;
        }
        //loglik += R::dpois(dist, rate, 1);
      }
    }
  }
  
  //Rcout << "Zero time counter = " << zero_time_counter << ", import length = " << import_idx.length() << " gen source = " 
  //      << gen_source << std::endl;
  return loglik;
}



double CalculateTransmissionLikelihood(List data, NumericVector parameters)
{
  double loglik = 0.0;
  bool print_details = false;
  int num_patients = data["num_patients"];
  IntegerVector t_a = data["t_a"];
  IntegerVector t_c = data["t_c"];
  IntegerVector t_d = data["t_d"];
  IntegerVector source = data["source"];
  double beta_p = parameters[2];
  double beta_h = 0.0;
  if(num_patients != t_c.length())
  {
    beta_h = parameters[3];
  }
  
  for(int i = 0; i<num_patients; i++)
  {
    if(t_c[i] == -1)
    {
      // The individual avoided infection the whole time
      for(int t = t_a[i]; t<=t_d[i]; t++)
      {
        loglik -= (beta_p*col_pop(data,t) + beta_h*hcw_col_pop(data,t));
      }
    }
    else if(source[i] == -1)
    {
      // The individual is an importation, and so no contribution to the likelihood
    }
    else if(t_c[i] == t_a[i])
    {
      // The individual was infected on their first day
      double total_infectious_pressure = beta_p*col_pop(data,t_c[i]) + beta_h*hcw_col_pop(data,t_c[i]);
      loglik += log(1-exp(-1*total_infectious_pressure)) - log(total_infectious_pressure);
      //loglik += log(1-exp(-(beta_p*col_pop(data,t_a[i]) + beta_h*hcw_col_pop(data,t_a[i]))));// - log(col_pop(data,t_c[i])+hcw_col_pop(data,t_c[i]));
      if(source[i] < num_patients)
      {
        //loglik += -1*log(col_pop(data,t_c[i])) + log(beta_p*col_pop(data,t_c[i])) - log(beta_p*col_pop(data,t_c[i]) + beta_h*hcw_col_pop(data,t_c[i]));
        loglik += log(beta_p);
      }
      else
      {
        //loglik += -1*log(hcw_col_pop(data,t_c[i])) + log(beta_h*hcw_col_pop(data,t_c[i])) - log(beta_p*col_pop(data,t_c[i]) + beta_h*hcw_col_pop(data,t_c[i]));
        loglik += log(beta_h);
      }
      
    }
    else
    {
      // The individual was infected after the first day
      for(int t = t_a[i]; t<t_c[i]; t++)
      {
        loglik -= (beta_p*col_pop(data,t) + beta_h*hcw_col_pop(data,t));
        //Rcout << "Day = " << t << ", loglik = " << loglik << std::endl;
      }
      //loglik += log(1-exp(-(beta_p*col_pop(data,t_c[i]) + beta_h*hcw_col_pop(data,t_c[i]))));// - log(col_pop(data,t_c[i])+hcw_col_pop(data,t_c[i]));
      double total_infectious_pressure = beta_p*col_pop(data,t_c[i]) + beta_h*hcw_col_pop(data,t_c[i]);
      loglik += log(1-exp(-1*total_infectious_pressure)) - log(total_infectious_pressure);
      if(source[i] < num_patients)
      {
        //loglik += -1*log(col_pop(data,t_c[i])) + log(beta_p*col_pop(data,t_c[i])) - log(beta_p*col_pop(data,t_c[i]) + beta_h*hcw_col_pop(data,t_c[i]));
        loglik += log(beta_p);
      }
      else
      {
        //loglik += -1*log(hcw_col_pop(data,t_c[i])) + log(beta_h*hcw_col_pop(data,t_c[i])) - log(beta_p*col_pop(data,t_c[i]) + beta_h*hcw_col_pop(data,t_c[i]));
        loglik += log(beta_h);
      }
    }
    if(print_details)
    {
      Rcout << "i = " << i << ", loglik = " << loglik << ", col pop = " << col_pop(data,t_c[i]) 
            << ", hcw col pop " << hcw_col_pop(data,t_c[i]) << ", source = " << source[i] << std::endl;
    }
    
    
  }
  return loglik;
}



// [[Rcpp::export]]
double CalculateLogLikelihood(List data, NumericVector parameters)
{
  bool print_details = false;
  double screening_likelihood = CalculateScreeningLikelihood(data, parameters);
  double importation_likelihood = CalculateImportationLikelihood(data, parameters);
  double transmission_likelihood = CalculateTransmissionLikelihood(data, parameters);
  double loglik = screening_likelihood + importation_likelihood + transmission_likelihood;
  double genetic_likelihood = 0.0;
  if(data.containsElementNamed("genetic_ids"))
  {
    genetic_likelihood = CalculateGeneticLikelihood(data, parameters);
    loglik += genetic_likelihood;
  }
  if(print_details)
  {
    Rcout << "Screen loglik = " << screening_likelihood << ", import loglik = " << importation_likelihood << ", transmission loglik = " 
          << transmission_likelihood << ", genetic loglik = " << genetic_likelihood << std::endl;
  }
  return loglik;
}

// [[Rcpp::export]]
double CalculateLogLikelihood_R(NumericVector parameters, IntegerVector t_a, IntegerVector t_c, IntegerVector t_d, IntegerVector source, 
                                IntegerVector hcw_ind, IntegerMatrix screening_matrix)
{
  
  List data = InitialiseData(t_a, t_c, t_d, source, hcw_ind, screening_matrix);
  bool print_details = true;
  double screening_likelihood = CalculateScreeningLikelihood(data, parameters);
  double importation_likelihood = CalculateImportationLikelihood(data, parameters);
  double transmission_likelihood = CalculateTransmissionLikelihood(data, parameters);
  double loglik = screening_likelihood + importation_likelihood + transmission_likelihood;
  if(print_details)
  {
    Rcout << "Screen loglik = " << screening_likelihood << ", import loglik = " << importation_likelihood << ", transmission loglik = " 
          << transmission_likelihood << std::endl;
  }
  return loglik;
}


// [[Rcpp::export]]
IntegerVector ReturnGenSourceVector(List data)
{
  //Rcout << "Begin calculating gen source\n";
  //Rcout << data << std::endl;
  //Rcout << "Before gen source" << std::endl;
  IntegerVector genetic_ids = data["genetic_ids"];
  IntegerVector source_vector = data["source"];
  IntegerVector sample_times = data["sample_times"];
  IntegerVector variant_numbers = data["variant_numbers"];
  IntegerVector t_c = data["t_c"];
  IntegerVector gen_source_vector;
  for(int i = 0; i<genetic_ids.length(); i++)
  {
    int source = ReturnSourceSequence(i, data);
    gen_source_vector.push_back(source);
  }
  //Rcout << "End calculating gen source\n";
  //Rcout << "After gen source" << std::endl;
  return gen_source_vector;
}

// [[Rcpp::export]]
int ReturnSourceSequence(int sequence_loc, List data)
{
  bool print_debug = false;
  
  IntegerVector genetic_ids = data["genetic_ids"];
  IntegerVector source_vector = data["source"];
  IntegerVector sample_times = data["sample_times"];
  IntegerVector variant_numbers = data["variant_numbers"];
  IntegerVector t_c = data["t_c"];
  
  if(print_debug) Rcout << "----- Sequence loc = " << sequence_loc << "\n";
  
  // First check if there are any other sequences at the t
  
  // first check if the sequence is a duplicate
  IntegerVector duplicated_sequences = ReturnDuplicateSequences(sequence_loc, data);
  if(print_debug) Rcout << "Duplicated sequences = " << duplicated_sequences << "\n";
  if(duplicated_sequences.length() > 0)
  {
    // there exists duplicate sequences, enforce each one points to the previous one, unless in position zero, then it must look backwards
    int sequence_idx = WhichVec(sequence_loc, duplicated_sequences)[0];
    if(sequence_idx != 0)
    {
      // return the previous one
      if(print_debug)
      {
        Rcout << "Return duplicated sequence idx = " << sequence_idx-1 << std::endl;
      }
      return duplicated_sequences[sequence_idx-1];
    }
  }
  
  // now look back to find the most recent swab or infection
  int current_variant = variant_numbers[sequence_loc];
  IntegerVector variant_idx = WhichVec(current_variant, variant_numbers);
  int target = genetic_ids[sequence_loc];
  int sequence_time = sample_times[sequence_loc];
  bool sequence_found = false;
  int source_sequence = -1;
  if(print_debug) Rcout << "Target = " << target << ", sequence time = " << sequence_time << "\n";
  
  // while the sequence is not found, look at the targets swabs and infection swabs at the time of infection until a sequence is found
  
  // cowboy workaround, this needs to be fixed asap
  //int time_adjust_counter = 0;
  while(!sequence_found)
  {
    // check the swabs
    IntegerVector target_swabs = WhichVec(target, genetic_ids);
    IntegerVector target_swabs_same_variant = intersect(target_swabs, variant_idx);
    IntegerVector target_swab_times_same_variant = as<IntegerVector>(sample_times[target_swabs_same_variant]);
    if(print_debug) Rcout << "Target = " << target << ", swabs = " << target_swabs_same_variant << ", times = " << target_swab_times_same_variant << "\n";
    int most_recent_time = -1;
    int most_recent_source = -1;
    for(int i = 0; i<target_swabs_same_variant.length(); i++)
    {
      int swab = target_swabs_same_variant[i];
      int swab_time = target_swab_times_same_variant[i];
      IntegerVector swab_loc_in_duplicates = WhichVec(swab, duplicated_sequences);
      if(swab_loc_in_duplicates.length() == 0)
      {
        if(swab_time > most_recent_time && swab_time < sequence_time)
        {
          most_recent_time = swab_time;
          most_recent_source = swab;
          // TESTING TO REDUCE TIME
          if(most_recent_time == sequence_time - 1)
          {
            sequence_found = true;
            break;
          }
        }
      }
      
    }
    
    if(sequence_found)
    {
      source_sequence = most_recent_source;
      break;
    }
    
    // now check the infections at the same time
    IntegerVector target_infections = WhichVec(target, source_vector);
    if(print_debug) Rcout << "Target infections = " << target_infections << "\n";
    int most_recent_inf_time = -1;
    int most_recent_inf_source = -1;
    for(int i = 0; i<target_infections.length(); i++)
    {
      int current_target = target_infections[i];
      int target_inf_time = t_c[current_target];
      // check if they were infected at the sequence time, if so then check for swabs
      IntegerVector inf_target_swabs = WhichVec(current_target, genetic_ids);
      if(print_debug) Rcout << "Infection target swabs = " << inf_target_swabs << "\n";
      IntegerVector inf_target_swabs_same_variant = intersect(inf_target_swabs, variant_idx);
      IntegerVector inf_target_swab_times_same_variant = as<IntegerVector>(sample_times[inf_target_swabs_same_variant]);
      IntegerVector loc_on_inf = WhichVec(target_inf_time, inf_target_swab_times_same_variant);
      if(print_debug)
      {
        Rcout << "Current infection target = " << current_target << ", infection time = " << target_inf_time << "\n"
              << "Infection target swabs times = " << inf_target_swab_times_same_variant << ", loc on infection = " << loc_on_inf << "\n";
      }
      if(loc_on_inf.length()>0)
      {
        int first_loc_on_inf = loc_on_inf[0];
        if(target_inf_time > most_recent_inf_time && target_inf_time < sequence_time && inf_target_swabs_same_variant[first_loc_on_inf] != sequence_loc)
        {
          IntegerVector loc_in_dup_seq = WhichVec(inf_target_swabs_same_variant[first_loc_on_inf],duplicated_sequences);
          if(loc_in_dup_seq.length()==0)
          {
            most_recent_inf_source = inf_target_swabs_same_variant[first_loc_on_inf];
            most_recent_inf_time = target_inf_time;
            // TESTING TO REDUCE TIME
            if(most_recent_inf_time == sequence_time - 1)
            {
              break;
            }
          }
          
          
        }
      }
    }
    if(print_debug) Rcout << "Most recent swab time = " << most_recent_time << ", most recent infection swab time = " << most_recent_inf_time << "\n";
    
    
    // check which was more recent, the swab or the swab of infections
    if(most_recent_time > most_recent_inf_time)
    {
      source_sequence = most_recent_source;
    }
    else
    {
      source_sequence = most_recent_inf_source;
    }
    
    IntegerVector source_loc_in_duplicated = WhichVec(source_sequence, duplicated_sequences);
    //if(source_loc_in_duplicated.length()>0) source_sequence = -1;
    
    if(source_sequence != -1)
    {
      sequence_found = true;
    }
    else
    {
      // couldnt find a sequence in the current target, look at the source
      sequence_time = t_c[target] + 1; // + 1 to account for at the same time
      target = source_vector[target];
      
      
      if(target==-1)
      {
        return -1;
      }
    }
  }
  
  if(print_debug) Rcout << "Try to see if the sequence is a duplicate" << std::endl;
  
  // now a sequence has been found, check if there are duplicates
  IntegerVector source_duplicates = ReturnDuplicateSequences(source_sequence, data);
  if(source_duplicates.length() > 0)
  {
    // there are duplicates of the source sequence, target the last sequence
    if(print_debug) Rcout << "Source duplicates = " << source_duplicates << "\n";
    return source_duplicates[source_duplicates.length()-1];
  }
  if(print_debug) Rcout << "Returning source sequence = " << source_sequence << std::endl;
  return source_sequence;
}

IntegerVector ReturnDuplicateSequences(int sequence_loc, List data)
{
  // 
  bool print_details = false;
  
  IntegerVector genetic_ids = data["genetic_ids"];
  IntegerVector source_vector = data["source"];
  IntegerVector sample_times = data["sample_times"];
  IntegerVector variant_numbers = data["variant_numbers"];
  IntegerVector t_c = data["t_c"];
  
  // create a vector to store all the duplicates
  IntegerVector duplicated_sequences;
  
  // first check if the sequence is at the time of infection
  int patient_id = genetic_ids[sequence_loc];
  int seq_sample_time = sample_times[sequence_loc];
  int current_variant = variant_numbers[sequence_loc];
  IntegerVector variant_loc = WhichVec(current_variant, variant_numbers);
  //int inf_time = t_c[patient_id];
  
  // check if there are 
  // (i) any other swabs at the time from the same patient
  // (ii) anyone who was infected by the target at the same time and also has a swab at that time (will only happen for imports)
  // (iii) any swabs from the source of infector
  // (iv) any swabs from other infections from the source
  
  if(print_details)
  {
    Rcout << "Variant locs = " << variant_loc << std::endl;
  }
  
  
  // check other swabs (i)
  IntegerVector patient_swabs = WhichVec(patient_id, genetic_ids);
  IntegerVector patient_variant_swabs = intersect(patient_swabs, variant_loc);
  IntegerVector patient_times = as<IntegerVector>(sample_times[patient_variant_swabs]);
  IntegerVector swabs_loc_at_time = WhichVec(seq_sample_time, patient_times);
  if(print_details)
  {
    Rcout << "Patient swabs = " << patient_swabs << std::endl;
    Rcout << "Patient variant swabs = " << patient_variant_swabs << std::endl;
    Rcout << "Patient times = " << patient_times << std::endl;
    Rcout << "Swabs loc at time = " << swabs_loc_at_time << std::endl;
  }
  
  for(int i = 0; i<swabs_loc_at_time.length(); i++)
  {
    int current_seq_idx = patient_variant_swabs[swabs_loc_at_time[i]];
    duplicated_sequences.push_back(current_seq_idx);
  }
  
  // check anyone who was infected by the target at the same time
  IntegerVector patient_infections = WhichVec(patient_id, source_vector);
  for(int i = 0; i<patient_infections.length(); i++)
  {
    int current_patient_infection = patient_infections[i];
    int patient_infection_col_time = t_c[current_patient_infection];
    if(patient_infection_col_time == seq_sample_time)
    {
      // an infection from the patient is at the same time the sequence was taken, now check if they have a swab
      IntegerVector target_swabs = WhichVec(current_patient_infection, genetic_ids);
      IntegerVector target_variant_swabs = intersect(target_swabs, variant_loc);
      IntegerVector target_times = as<IntegerVector>(sample_times[target_variant_swabs]);
      IntegerVector target_swabs_loc_at_time = WhichVec(seq_sample_time, target_times);
      for(int j = 0; j<target_swabs_loc_at_time.length(); j++)
      {
        int current_seq_idx = target_variant_swabs[target_swabs_loc_at_time[j]];
        duplicated_sequences.push_back(current_seq_idx);
      }
    }
  }
  
  if(print_details)
  {
    Rcout << "After swabs, dup sequences = " << duplicated_sequences << std::endl;
  }
  
  // check swabs from the source of infection, only if t_c = sequence time
  int patient_col_time = t_c[patient_id];
  
  if(patient_col_time == seq_sample_time)
  {
    if(print_details) Rcout << "Checking the source of patient " << patient_id << std::endl;
    int patient_source = source_vector[patient_id];
    if(print_details)
    {
      Rcout << "patient source = " << patient_source << std::endl;
    }
    if(patient_source != -1)
    {
      // if they are not an importation
      IntegerVector target_swabs = WhichVec(patient_source, genetic_ids);
      IntegerVector target_variant_swabs = intersect(target_swabs, variant_loc);
      IntegerVector target_times = as<IntegerVector>(sample_times[target_variant_swabs]);
      IntegerVector target_swabs_loc_at_time = WhichVec(seq_sample_time, target_times);
      if(print_details)
      {
        Rcout << "Target swabs = " << target_swabs << std::endl;
        Rcout << "Target times = " << target_times << std::endl;
        Rcout << "target_swabs_loc_at_time = " << target_swabs_loc_at_time << std::endl;
      }
      for(int j = 0; j<target_swabs_loc_at_time.length(); j++)
      {
        int current_seq_idx = target_variant_swabs[target_swabs_loc_at_time[j]];
        duplicated_sequences.push_back(current_seq_idx);
      }
    }
    
    if(print_details)
    {
      Rcout << "After swabs from source, dup sequences = " << duplicated_sequences << std::endl;
    }
    
    // now check other infections from the source
    if(patient_source != -1)
    {
      IntegerVector source_infections = WhichVec(patient_source, source_vector);
      for(int i = 0; i<source_infections.length(); i++)
      {
        int current_infection_id = source_infections[i];
        if(current_infection_id != patient_id)
        {
          int current_infection_time = t_c[current_infection_id];
          if(current_infection_time == seq_sample_time)
          {
            // there is another infection at the same time, check their swabs
            IntegerVector target_swabs = WhichVec(current_infection_id, genetic_ids);
            IntegerVector target_variant_swabs = intersect(target_swabs, variant_loc);
            IntegerVector target_times = as<IntegerVector>(sample_times[target_variant_swabs]);
            IntegerVector target_swabs_loc_at_time = WhichVec(seq_sample_time, target_times);
            for(int j = 0; j<target_swabs_loc_at_time.length(); j++)
            {
              int current_seq_idx = target_variant_swabs[target_swabs_loc_at_time[j]];
              duplicated_sequences.push_back(current_seq_idx);
            }
          }
        }
      }
    }
    
    
    if(print_details)
    {
      Rcout << "After swabs from source infections, dup sequences = " << duplicated_sequences << std::endl;
    }
  }
  
  
  
  
  
  if(duplicated_sequences.length()>1)
  {
    if(print_details) Rcout << "Returning duplicate sequences = " << sort_unique(duplicated_sequences) << std::endl;
    return sort_unique(duplicated_sequences);
  }
  else
  {
    IntegerVector x;
    if(print_details) Rcout << "Returning duplicate sequences = " << x << std::endl;
    return x;
  }
  
  
}

// [[Rcpp::export]]
IntegerVector Rcpp_sort(IntegerVector x, IntegerVector y) {
  // Order the elements of x by sorting y
  // First create a vector of indices
  IntegerVector idx = seq_along(x) - 1;
  // Then sort that vector by the values of y
  std::sort(idx.begin(), idx.end(), [&](int i, int j){return y[i] < y[j];});
  // And return x in that order
  return x[idx];
}

int CalculateImportationSum(List data)
{
  int num_patients = data["num_patients"];
  IntegerVector source = data["source"];
  int import_sum = 0;
  for(int i = 0; i<num_patients; i++)
  {
    if(source[i] == -1)
    {
      import_sum++;
    }
  }
  return import_sum;
}

List CalculateImportDistance(List data, NumericVector parameters) {
  IntegerVector gen_source = ReturnGenSourceVector(data);
  IntegerVector source = data["source"];
  IntegerVector variant_numbers = data["variant_numbers"];
  IntegerMatrix genetic_matrix = data["genetic_matrix"];
  int distance_sum = 0;
  int variant_distance_sum = 0;
  double import_counter = 0;
  double variant_import_counter = 0;
  
  IntegerVector import_idx = WhichVec(-1, gen_source);
  //Rcout << "Import idx = " << import_idx << std::endl;
  for(int j = 1; j<import_idx.length(); j++)
  {
    for(int i = 0; i<j; i++)
    {
      int ii = import_idx[i];
      int jj = import_idx[j];
      if(variant_numbers[ii] == variant_numbers[jj]) 
      {
        //Rcout << "genetic_matrix(ii,jj) = " << genetic_matrix(ii,jj) << std::endl;
        distance_sum += genetic_matrix(ii,jj);
        import_counter++;
      }
      else
      {
        variant_distance_sum += genetic_matrix(ii,jj);
        variant_import_counter++;
      }
    }
  }
  
  List out = List::create(Named("distance_sum") = distance_sum,
                          Named("import_counter") = import_counter,
                          Named("variant_distance_sum") = variant_distance_sum,
                          Named("variant_import_counter") = variant_import_counter);
  return out;
}


// Vq is the set of patients who tested positive during their stay
IntegerVector CalculateVq(List data)
{
  IntegerVector Vq;
  IntegerMatrix screening_matrix = data["screening_matrix"];
  int num_patients = data["num_patients"];
  for(int i = 0; i<num_patients; i++)
  {
    IntegerVector current_row = screening_matrix(i,_);
    IntegerVector positive_locs = WhichVec(1, current_row);
    if(positive_locs.length()>0)
    {
      Vq.push_back(i);
    }
  }
  Rcout << "Calculated Vq" << std::endl;
  return Vq;
}

IntegerVector CalculateVs(List data)
{
  IntegerVector Vs;
  IntegerMatrix screening_matrix = data["screening_matrix"];
  int num_patients = data["num_patients"];
  for(int i = 0; i<num_patients; i++)
  {
    IntegerVector current_row = screening_matrix(i,_);
    IntegerVector positive_locs = WhichVec(1, current_row);
    if(positive_locs.length()==0)
    {
      Vs.push_back(i);
    }
  }
  Rcout << "Calculated Vs" << std::endl;
  return Vs;
}


void UpdateTransmissionRate_P(List data, NumericVector &parameters, double &loglik_cur, double proposal_variance, double prior, int &nacc)
{
  bool print_details = false;
  NumericVector parameters_can = clone(parameters);
  parameters_can[2] = R::rnorm(parameters[2], proposal_variance);
  if(parameters_can[2] > 0)
  {
    double logpi_cur = loglik_cur + R::dexp(parameters[2], 1/prior, 1);
    double loglik_can = CalculateLogLikelihood(data, parameters_can);
    double logpi_can = loglik_can + R::dexp(parameters_can[2], 1/prior, 1);
    double U = R::runif(0.0,1.0);
    if(print_details)
    {
      Rcout << "Beta cur = " << parameters[2] << ", beta can = " << parameters_can[2]
            << "loglik cur = " << loglik_cur << ", loglik can = " << loglik_can;
    }
    if(log(U) < logpi_can - logpi_cur)
    {
      if(print_details)
      {
        Rcout << " - accepted." << std::endl;
      }
      parameters = parameters_can;
      loglik_cur = loglik_can;
      nacc++;
    }
    else
    {
      if(print_details)
      {
        Rcout << " - rejected." << std::endl;
      }
    }
  }
}

void UpdateTransmissionRate_H(List data, NumericVector &parameters, double &loglik_cur, double proposal_variance, double prior, int &nacc)
{
  bool print_details = false;
  NumericVector parameters_can = clone(parameters);
  parameters_can[3] = R::rnorm(parameters[3], proposal_variance);
  if(parameters_can[3] > 0)
  {
    double logpi_cur = loglik_cur + R::dexp(parameters[3], 1/prior, 1);
    double loglik_can = CalculateLogLikelihood(data, parameters_can);
    double logpi_can = loglik_can + R::dexp(parameters_can[3], 1/prior, 1);
    double U = R::runif(0.0,1.0);
    if(print_details)
    {
      Rcout << "Beta cur = " << parameters[3] << ", beta can = " << parameters_can[3]
            << "loglik cur = " << loglik_cur << ", loglik can = " << loglik_can;
    }
    if(log(U) < logpi_can - logpi_cur)
    {
      if(print_details)
      {
        Rcout << " - accepted." << std::endl;
      }
      parameters = parameters_can;
      loglik_cur = loglik_can;
      nacc++;
    }
    else
    {
      if(print_details)
      {
        Rcout << " - rejected." << std::endl;
      }
    }
  }
}

void UpdateMutationRate(List data, NumericVector &parameters, double &loglik_cur, double proposal_variance, double prior, int &nacc)
{
  NumericVector parameters_can = clone(parameters);
  parameters_can[4] = R::rnorm(parameters[4], proposal_variance);
  if(parameters_can[4] > 0)
  {
    double logpi_cur = loglik_cur + R::dexp(parameters[4], 1/prior, 1);
    double loglik_can = CalculateLogLikelihood(data, parameters_can);
    double logpi_can = loglik_can + R::dexp(parameters_can[4], 1/prior, 1);
    double U = R::runif(0.0,1.0);
    if(log(U) < logpi_can - logpi_cur)
    {
      parameters = parameters_can;
      loglik_cur = loglik_can;
      nacc++;
    }
  }
}

int ReturnLastDay(List data, int target)
{
  // A targets last day is either, the day of discharge, the day of the first positive test result, or the day of the first offspring birth
  // whatever comes first
  IntegerVector t_a = data["t_c"];
  IntegerVector t_c = data["t_c"];
  IntegerVector t_d = data["t_d"];
  IntegerVector source = data["source"];
  IntegerMatrix screening_matrix = data["screening_matrix"];
  int last_day = t_d[target];
  //Rcout << "Target = " << target << "t_a = " << t_a[target] << ", offspring coltime = ";
  IntegerVector offspring = WhichVec(target, source);
  for(int i = 0; i<offspring.length(); i++)
  {
    int offspring_coltime = t_c[offspring[i]];
    if(offspring_coltime <= last_day)
    {
      last_day = offspring_coltime - 1;// - 1;  // minus one because you can only infect someone on the day before
    }
  }
  
  IntegerVector screening_row = screening_matrix(target,_);
  IntegerVector positive_days = WhichVec(1,screening_row);
  //Rcout << "positive days = " << positive_days << std::endl;
  if(positive_days.length()>0)
  {
    int first_swab_day = positive_days[0];
    if(first_swab_day < last_day)
    {
      last_day = first_swab_day;
    }
  }
  //Rcout << ", last day = " << last_day << std::endl;
  
  return last_day;
}

IntegerVector ReturnNadd0(List data, IntegerVector Va)
{
  IntegerVector source = data["source"];
  IntegerVector nadd0;
  for(int i = 0; i<Va.length(); i++)
  {
    IntegerVector infections = WhichVec(Va[i], source);
    if(infections.length()==0)
    {
      nadd0.push_back(Va[i]);
    }
  }
  return nadd0;
}

IntegerVector ReturnPossibleInfectors(List data, int target)
{
  IntegerVector t_c = data["t_c"];
  IntegerVector t_d = data["t_d"];
  IntegerVector source = data["source"];
  IntegerVector possible_infectors;
  
  int num_patients = data["num_patients"];
  int day_of_colonisation = t_c[target];
  for(int i = 0; i<num_patients; i++)
  {
    if(t_c[i] != -1 && t_c[i] <= day_of_colonisation && day_of_colonisation <= t_d[i] && i != target)
    {
      if(t_c[i] == day_of_colonisation && source[i] == -1)
      {
        possible_infectors.push_back(i);
      }
      else if(t_c[i] < day_of_colonisation)
      {
        possible_infectors.push_back(i);
      }
    }
  }
  
  // Now check HCWs
  for(int i = num_patients; i<t_c.length(); i++)
  {
    if(t_c[i] != -1 && t_c[i] < day_of_colonisation && day_of_colonisation <= t_d[i] && i != target)
    {
      possible_infectors.push_back(i);
    }
  }
  
  
  return possible_infectors;
}

IntegerVector ReturnPossibleInfectors_patient(List data, int target)
{
  IntegerVector t_c = data["t_c"];
  IntegerVector t_d = data["t_d"];
  IntegerVector source = data["source"];
  IntegerVector possible_infectors;
  
  int num_patients = data["num_patients"];
  int day_of_colonisation = t_c[target];
  for(int i = 0; i<num_patients; i++)
  {
    if(t_c[i] != -1 && t_c[i] <= day_of_colonisation && day_of_colonisation <= t_d[i] && i != target)
    {
      if(t_c[i] == day_of_colonisation && source[i] == -1)
      {
        possible_infectors.push_back(i);
      }
      else if(t_c[i] < day_of_colonisation)
      {
        possible_infectors.push_back(i);
      }
    }
  }
  
  return possible_infectors;
}

IntegerVector ReturnPossibleInfectors_hcw(List data, int target)
{
  IntegerVector t_c = data["t_c"];
  IntegerVector t_d = data["t_d"];
  IntegerVector source = data["source"];
  IntegerVector possible_infectors;
  
  int num_patients = data["num_patients"];
  int day_of_colonisation = t_c[target];
  // Now check HCWs
  for(int i = num_patients; i<t_c.length(); i++)
  {
    if(t_c[i] != -1 && t_c[i] < day_of_colonisation && day_of_colonisation <= t_d[i] && i != target)
    {
      possible_infectors.push_back(i);
    }
  }
  
  return possible_infectors;
}

List MoveColonisationTimeGenetic(List &data, NumericVector parameters, double &loglik_cur, double w, IntegerVector Vq, IntegerVector Va, int &nacc_move, std::ofstream &myfile)
{
  bool print_details = false;
  bool write_details_to_file = true;
  
  IntegerVector t_a = data["t_a"];
  IntegerVector t_d = data["t_d"];
  IntegerVector t_c = data["t_c"];
  IntegerVector source = data["source"];
  
  List data_can = clone(data);
  IntegerVector t_c_can = data_can["t_c"];
  IntegerVector source_can = data_can["source"];
  
  IntegerVector colonised_population = union_(Vq,Va);
  int target = SampleVector(colonised_population);
  int last_day = ReturnLastDay(data, target);
  
  double log_prop_ratio = 0.0;
  double U = R::runif(0.0,1.0);
  if(U < w)
  {
    // Propose the target is an importation
    t_c_can[target] = t_a[target];
    source_can[target] = -1;
    if(source[target]==-1)
    {
      // Importation -> Importation
      // Do nothing, proposal ratio = 1
    }
    else
    {
      // Acquisition -> Importation
      log_prop_ratio = log(1-w) - log(w) - log(last_day-t_a[target]+1) - log(total_col_pop(data, t_c[target]));
    }
  }
  else
  {
    // Propose the target is an acquisition
    if(last_day < t_a[target]) return data; // The target cannot be an acquisition
    t_c_can[target] = SampleVector(seq(t_a[target],last_day));
    IntegerVector possible_infectors = ReturnPossibleInfectors(data_can, target);
    //Rcout << "Possible infectors " << possible_infectors << std::endl;
    if(possible_infectors.length()==0) return data;
    source_can[target] = SampleVector(possible_infectors);
    if(source[target] != -1)
    {
      // Acquisition -> Acquisition
      log_prop_ratio = log(total_col_pop(data_can, t_c_can[target])) - log(total_col_pop(data, t_c[target]));
    }
    else
    {
      // Importation -> Acquisition
      log_prop_ratio = log(w) + log(last_day - t_a[target] + 1) - log(1-w) 
                        + log(total_col_pop(data_can, t_c_can[target]));
    }
  }
  
  //if(t_c[target] == t_c_can[target] && source[target] == source_can[target]) return;
  
  bool impute_distances = true;
  
  if(data.containsElementNamed("genetic_ids") && impute_distances)
  {
    //data_can = UpdateImputedNodesMove3(data, data_can, parameters, target);
    data_can = UpdateImputedNodesMoveAll(data, data_can, parameters);
    double genetic_contribution = data_can["genetic_contribution"];
    if(genetic_contribution < -500000 || genetic_contribution > 5000000)
    {
      Rcout << "Genetic contribution = " << genetic_contribution << std::endl;
      stop("problematic genetic contribution");
    }
    log_prop_ratio += genetic_contribution;
  }
  
  
  double loglik_can = CalculateLogLikelihood(data_can, parameters);
  U = R::runif(0.0,1.0);
  if(print_details)
  {
    Rcout << "MOVE ID = " << target << " from time t = " << t_c[target] << " to t = " << t_c_can[target]
          << " from source " << source[target] << " to source " << source_can[target] << std::endl;
    Rcout << "Loglik can = " << loglik_can << ", loglik cur = " << loglik_cur << ", log prop ratio = " << log_prop_ratio;
  }
  
  if(write_details_to_file)
  {
    myfile << "MOVE ID = " << target << " from time t = " << t_c[target] << " to t = " << t_c_can[target]
           << " from source " << source[target] << " to source " << source_can[target] << std::endl;
    myfile << "Loglik can = " << loglik_can << ", loglik cur = " << loglik_cur << ", log prop ratio = " << log_prop_ratio << std::endl;
    
    double transmission_ll_cur = CalculateTransmissionLikelihood(data, parameters);
    double transmission_ll_can = CalculateTransmissionLikelihood(data_can, parameters);
    
    double screening_ll_cur = CalculateScreeningLikelihood(data, parameters);
    double screening_ll_can = CalculateScreeningLikelihood(data_can, parameters);
    
    double import_ll_cur = CalculateImportationLikelihood(data, parameters);
    double import_ll_can = CalculateImportationLikelihood(data_can, parameters);
    
    double genetic_ll_cur = CalculateGeneticLikelihood(data, parameters);
    double genetic_ll_can = CalculateGeneticLikelihood(data_can, parameters);
    
    myfile << "screen ll cur = " << screening_ll_cur << ", import ll cur = " << import_ll_cur << ", transmission ll cur = "
           << transmission_ll_cur << ", genetic ll cur = " << genetic_ll_cur << std::endl;
    myfile << "screen ll can = " << screening_ll_can << ", import ll can = " << import_ll_can << ", transmission ll can = "
           << transmission_ll_can << ", genetic ll can = " << genetic_ll_can << std::endl;
    
    IntegerVector gen_source_cur = ReturnGenSourceVector(data);
    IntegerVector gen_source_can = ReturnGenSourceVector(data_can);
    
    myfile << "       Gen source cur = ";
    
    for(int i = 0; i<gen_source_cur.length(); i++)
    {
      int current_gen_source = gen_source_cur[i];
      if(current_gen_source == -1)
      {
        myfile << current_gen_source << " ";
      }
      else if(current_gen_source < 10)
      {
        myfile << " " << current_gen_source << " ";
      }
      else
      {
        myfile << current_gen_source << " ";
      }
    }
    myfile << std::endl;
    
    
    myfile << "       Gen source can = ";
    for(int i = 0; i<gen_source_can.length(); i++)
    {
      int current_gen_source = gen_source_can[i];
      if(current_gen_source == -1)
      {
        myfile << current_gen_source << " ";
      }
      else if(current_gen_source < 10)
      {
        myfile << " " << current_gen_source << " ";
      }
      else
      {
        myfile << current_gen_source << " ";
      }
    }
    myfile << std::endl;
    
    IntegerMatrix gen_matrix_cur = data["genetic_matrix"];
    IntegerMatrix gen_matrix_can = data_can["genetic_matrix"];
    
    myfile << "observed distance cur = ";
    for(int i = 0; i<gen_source_cur.length(); i++)
    {
      int current_gen_source = gen_source_cur[i];
      int distance = gen_matrix_cur(i, current_gen_source);
      if(current_gen_source != -1)
      {
        if(distance < 10)
        {
          myfile << " " << distance << " ";
        }
        else
        {
          myfile << distance << " ";
        }
      }
      else
      {
        myfile << "-1 "; 
      }
    }
    myfile << std::endl;
    
    myfile << "observed distance can = ";
    for(int i = 0; i<gen_source_can.length(); i++)
    {
      int current_gen_source = gen_source_can[i];
      int distance = gen_matrix_can(i, current_gen_source);
      if(current_gen_source != -1)
      {
        if(distance < 10)
        {
          myfile << " " << distance << " ";
        }
        else
        {
          myfile << distance << " ";
        }
      }
      else
      {
        myfile << "-1 "; 
      }
    }
    myfile << std::endl;
    
    IntegerVector sample_times_cur = data["sample_times"];
    IntegerVector sample_times_can = data_can["sample_times"];
    myfile << "     sample times cur = ";
    for(int i = 0; i<sample_times_cur.length(); i++)
    {
      int current_gen_source = gen_source_cur[i];
      if(current_gen_source != -1)
      {
        int time_diff = sample_times_cur[i] - sample_times_cur[current_gen_source];
        if(time_diff < 10)
        {
          myfile << " " << time_diff << " ";
        }
        else
        {
          myfile << time_diff << " ";
        }
        
      }
      else
      {
        myfile << "-1 ";
      }
    }
    myfile << std::endl;
    
    myfile << "     sample times can = ";
    for(int i = 0; i<sample_times_can.length(); i++)
    {
      int current_gen_source = gen_source_can[i];
      if(current_gen_source != -1)
      {
        int time_diff = sample_times_can[i] - sample_times_can[current_gen_source];
        if(time_diff < 10)
        {
          myfile << " " << time_diff << " ";
        }
        else
        {
          myfile << time_diff << " ";
        }
      }
      else
      {
        myfile << "-1 ";
      }
    }
    myfile << std::endl;
    
    myfile << "- Current Transmission tree - " << std::endl;
    myfile << t_a << std::endl;
    myfile << t_c << std::endl;
    myfile << t_d << std::endl;
    myfile << source << std::endl;
    
    IntegerVector genetic_ids = data["genetic_ids"];
    IntegerVector sample_times = data["sample_times"];
    myfile << genetic_ids << std::endl;
    myfile << sample_times << std::endl;
  }
  
  
  if(log(U) < loglik_can - loglik_cur + log_prop_ratio)
  {
    //data = data_can;
    loglik_cur = loglik_can;
    nacc_move++;
    if(print_details) Rcout << " - accepted." << std::endl;
    if(write_details_to_file) myfile << " - accepted." << std::endl << std::endl;
    return data_can;
  }
  else
  {
    if(print_details) Rcout << " - rejected." << std::endl;
    if(write_details_to_file) myfile << " - rejected." << std::endl << std::endl;
    return data;
  }
  
}


List MoveMultipleColonisationTimeGenetic(List data, NumericVector parameters, double &loglik_cur, double w, IntegerVector Vq, IntegerVector Va, int &nacc_move, std::ofstream &myfile)
{
  bool print_details = false;
  bool write_details_to_file = true;
  
  IntegerVector t_a = data["t_a"];
  IntegerVector t_d = data["t_d"];
  IntegerVector t_c = data["t_c"];
  IntegerVector source = data["source"];
  
  List data_cur = clone(data);
  
  List data_can = clone(data);
  IntegerVector t_c_can = data_can["t_c"];
  IntegerVector source_can = data_can["source"];
  
  IntegerVector colonised_population = union_(Vq,Va);
  //IntegerVector targets = SampleMultipleVector(colonised_population, 3);
  IntegerVector targets = {62,59,56};
  double log_prop_ratio = 0.0;
  
  for(int i = 0; i<targets.length(); i++)
  {
    int target = targets[i];
    int last_day = ReturnLastDay(data_cur, target);
    double U = R::runif(0.0,1.0);
    if(U < w)
    {
      // Propose the target is an importation
      t_c_can[target] = t_a[target];
      source_can[target] = -1;
      if(source[target]==-1)
      {
        // Importation -> Importation
        // Do nothing, proposal ratio = 1
      }
      else
      {
        // Acquisition -> Importation
        log_prop_ratio += log(1-w) - log(w) - log(last_day-t_a[target]+1) - log(total_col_pop(data_cur, t_c[target]));
      }
    }
    else
    {
      // Propose the target is an acquisition
      if(last_day < t_a[target]) return data; // The target cannot be an acquisition
      t_c_can[target] = SampleVector(seq(t_a[target],last_day));
      IntegerVector possible_infectors = ReturnPossibleInfectors(data_can, target);
      if(target == 56 && t_c_can[target] >= 68)
      {
        Rcout << "ID 62, t_c = " << t_c_can[62] << ", source = " << source_can[62] << std::endl;
        Rcout << "ID 56, t_c = " << t_c_can[target] << ", possible infectors = " << possible_infectors << std::endl;
        Rcout << "Last day of 56 = " << last_day << std::endl;
      }
      //Rcout << "Possible infectors " << possible_infectors << std::endl;
      if(possible_infectors.length()==0) return data;
      source_can[target] = SampleVector(possible_infectors);
      if(source[target] != -1)
      {
        // Acquisition -> Acquisition
        log_prop_ratio += log(total_col_pop(data_can, t_c_can[target])) - log(total_col_pop(data_cur, t_c[target]));
      }
      else
      {
        // Importation -> Acquisition
        log_prop_ratio += log(w) + log(last_day - t_a[target] + 1) - log(1-w) + log(total_col_pop(data_can, t_c_can[target]));
      }
    }
    data_cur = data_can;
  }


  
  //if(t_c[target] == t_c_can[target] && source[target] == source_can[target]) return;
  
  bool impute_distances = true;
  
  if(data.containsElementNamed("genetic_ids") && impute_distances)
  {
    //data_can = UpdateImputedNodesMove3(data, data_can, parameters, target);
    data_can = UpdateImputedNodesMoveAll(data, data_can, parameters);
    double genetic_contribution = data_can["genetic_contribution"];
    if(genetic_contribution < -500000 || genetic_contribution > 5000000)
    {
      Rcout << "Genetic contribution = " << genetic_contribution << std::endl;
      stop("problematic genetic contribution");
    }
    log_prop_ratio += genetic_contribution;
  }
  
  
  double loglik_can = CalculateLogLikelihood(data_can, parameters);
  double U = R::runif(0.0,1.0);
  if(print_details)
  {
    Rcout << "MOVE TARGETS = " << targets << std::endl;
    Rcout << "Loglik can = " << loglik_can << ", loglik cur = " << loglik_cur << ", log prop ratio = " << log_prop_ratio;
  }
  
  if(write_details_to_file)
  {
    myfile << "MOVE TARGETS = " << targets << std::endl;
    myfile << "Multi Proposal" << std::endl;
    myfile << "Target = 56, t_c_can = " << t_c_can[56] << ", source_can = " << source_can[56] << std::endl;
    myfile << "Target = 59, t_c_can = " << t_c_can[59] << ", source_can = " << source_can[59] << std::endl;
    myfile << "Target = 62, t_c_can = " << t_c_can[62] << ", source_can = " << source_can[62] << std::endl;
    myfile << "Loglik can = " << loglik_can << ", loglik cur = " << loglik_cur << ", log prop ratio = " << log_prop_ratio << std::endl;
    
    myfile << "Target = 56, t_c = " << t_c[56] << ", source = " << source[56] << std::endl;
    myfile << "Target = 59, t_c = " << t_c[59] << ", source = " << source[59] << std::endl;
    myfile << "Target = 62, t_c = " << t_c[62] << ", source = " << source[62] << std::endl;
    
    double transmission_ll_cur = CalculateTransmissionLikelihood(data, parameters);
    double transmission_ll_can = CalculateTransmissionLikelihood(data_can, parameters);
    
    double screening_ll_cur = CalculateScreeningLikelihood(data, parameters);
    double screening_ll_can = CalculateScreeningLikelihood(data_can, parameters);
    
    double import_ll_cur = CalculateImportationLikelihood(data, parameters);
    double import_ll_can = CalculateImportationLikelihood(data_can, parameters);
    
    double genetic_ll_cur = CalculateGeneticLikelihood(data, parameters);
    double genetic_ll_can = CalculateGeneticLikelihood(data_can, parameters);
    
    myfile << "screen ll cur = " << screening_ll_cur << ", import ll cur = " << import_ll_cur << ", transmission ll cur = "
           << transmission_ll_cur << ", genetic ll cur = " << genetic_ll_cur << std::endl;
    myfile << "screen ll can = " << screening_ll_can << ", import ll can = " << import_ll_can << ", transmission ll can = "
           << transmission_ll_can << ", genetic ll can = " << genetic_ll_can << std::endl;
    
    IntegerVector gen_source_cur = ReturnGenSourceVector(data);
    IntegerVector gen_source_can = ReturnGenSourceVector(data_can);
    
    myfile << "Gen source cur = " << gen_source_cur << std::endl;
    myfile << "Gen source can = " << gen_source_can << std::endl;
    
    IntegerMatrix gen_matrix_cur = data["genetic_matrix"];
    IntegerMatrix gen_matrix_can = data_can["genetic_matrix"];
    
    myfile << "observed distance cur = ";
    for(int i = 0; i<gen_source_cur.length(); i++)
    {
      int current_gen_source = gen_source_cur[i];
      if(current_gen_source != -1)
      {
        myfile << gen_matrix_cur(i, current_gen_source) << " ";
      }
      else
      {
        myfile << "-1 "; 
      }
      
    }
    myfile << std::endl;
    
    myfile << "observed distance can = ";
    for(int i = 0; i<gen_source_can.length(); i++)
    {
      int current_gen_source = gen_source_can[i];
      if(current_gen_source != -1)
      {
        myfile << gen_matrix_can(i, current_gen_source) << " ";
      }
      else
      {
        myfile << "-1 ";
      }
    }
    myfile << std::endl;
    
    IntegerVector sample_times_cur = data["sample_times"];
    IntegerVector sample_times_can = data_can["sample_times"];
    myfile << "sample times cur = ";
    for(int i = 0; i<sample_times_cur.length(); i++)
    {
      int current_gen_source = gen_source_cur[i];
      if(current_gen_source != -1)
      {
        int time_diff = sample_times_cur[i] - sample_times_cur[current_gen_source];
        myfile << time_diff << " ";
      }
      else
      {
        myfile << "-1 ";
      }
    }
    myfile << std::endl;
    
    myfile << "sample times can = ";
    for(int i = 0; i<sample_times_can.length(); i++)
    {
      int current_gen_source = gen_source_can[i];
      if(current_gen_source != -1)
      {
        int time_diff = sample_times_can[i] - sample_times_can[current_gen_source];
        myfile << time_diff << " ";
      }
      else
      {
        myfile << "-1 ";
      }
    }
    myfile << std::endl;
  }
  
  
  if(log(U) < loglik_can - loglik_cur + log_prop_ratio)
  {
    //data = data_can;
    loglik_cur = loglik_can;
    nacc_move++;
    if(print_details) Rcout << " - accepted." << std::endl;
    if(write_details_to_file) myfile << " - accepted." << std::endl << std::endl;
    return data_can;
  }
  else
  {
    if(print_details) Rcout << " - rejected." << std::endl;
    if(write_details_to_file) myfile << " - rejected." << std::endl << std::endl;
    return data;
  }
  
}

// Update a source without changing the time of colonisation
List UpdateSource(List &data, NumericVector parameters, double &loglik_cur, double w, IntegerVector Vq, IntegerVector Va, int &nacc_move, std::ofstream &myfile)
{
  bool print_details = false;
  bool write_details_to_file = true;
  
  IntegerVector t_a = data["t_a"];
  IntegerVector t_d = data["t_d"];
  IntegerVector t_c = data["t_c"];
  IntegerVector source = data["source"];
  
  List data_can = clone(data);
  //IntegerVector t_c_can = data_can["t_c"];
  IntegerVector source_can = data_can["source"];
  
  IntegerVector colonised_population = union_(Vq,Va);
  int target = SampleVector(colonised_population);
  //int last_day = ReturnLastDay(data, target);
  
  double log_prop_ratio = 0.0;
  
  if(source[target] == -1) return data;
  
  IntegerVector possible_infectors = ReturnPossibleInfectors(data_can, target);
  if(possible_infectors.length()==0) return data;
  source_can[target] = SampleVector(possible_infectors);
  
  bool impute_distances = true;
  
  if(data.containsElementNamed("genetic_ids") && impute_distances)
  {
    //data_can = UpdateImputedNodesMove3(data, data_can, parameters, target);
    data_can = UpdateImputedNodesMoveAll(data, data_can, parameters);
    double genetic_contribution = data_can["genetic_contribution"];
    if(genetic_contribution < -500000 || genetic_contribution > 5000000)
    {
      Rcout << "Genetic contribution = " << genetic_contribution << std::endl;
      stop("problematic genetic contribution");
    }
    log_prop_ratio += genetic_contribution;
  }
  
  
  double loglik_can = CalculateLogLikelihood(data_can, parameters);
  double U = R::runif(0.0,1.0);
  if(print_details)
  {
    Rcout << "CHANGE SOURCE ID = " << target << " from source " << source[target] << " to source " << source_can[target] << std::endl;
    Rcout << "Loglik can = " << loglik_can << ", loglik cur = " << loglik_cur << ", log prop ratio = " << log_prop_ratio;
  }
  
  if(write_details_to_file)
  {
    myfile << "CHANGE SOURCE ID = " << target << " from source " << source[target] << " to source " << source_can[target] << std::endl;
    myfile << "Loglik can = " << loglik_can << ", loglik cur = " << loglik_cur << ", log prop ratio = " << log_prop_ratio << std::endl;
    
    double transmission_ll_cur = CalculateTransmissionLikelihood(data, parameters);
    double transmission_ll_can = CalculateTransmissionLikelihood(data_can, parameters);
    
    double screening_ll_cur = CalculateScreeningLikelihood(data, parameters);
    double screening_ll_can = CalculateScreeningLikelihood(data_can, parameters);
    
    double import_ll_cur = CalculateImportationLikelihood(data, parameters);
    double import_ll_can = CalculateImportationLikelihood(data_can, parameters);
    
    double genetic_ll_cur = CalculateGeneticLikelihood(data, parameters);
    double genetic_ll_can = CalculateGeneticLikelihood(data_can, parameters);
    
    myfile << "screen ll cur = " << screening_ll_cur << ", import ll cur = " << import_ll_cur << ", transmission ll cur = "
           << transmission_ll_cur << ", genetic ll cur = " << genetic_ll_cur << std::endl;
    myfile << "screen ll can = " << screening_ll_can << ", import ll can = " << import_ll_can << ", transmission ll can = "
           << transmission_ll_can << ", genetic ll can = " << genetic_ll_can << std::endl;
    
    IntegerVector gen_source_cur = ReturnGenSourceVector(data);
    IntegerVector gen_source_can = ReturnGenSourceVector(data_can);
    
    myfile << "Gen source cur = " << gen_source_cur << std::endl;
    myfile << "Gen source can = " << gen_source_can << std::endl;
    
    IntegerMatrix gen_matrix_cur = data["genetic_matrix"];
    IntegerMatrix gen_matrix_can = data_can["genetic_matrix"];
    
    myfile << "observed distance cur = ";
    for(int i = 0; i<gen_source_cur.length(); i++)
    {
      int current_gen_source = gen_source_cur[i];
      if(current_gen_source != -1)
      {
        myfile << gen_matrix_cur(i, current_gen_source) << " ";
      }
      else
      {
        myfile << "-1 "; 
      }
      
    }
    myfile << std::endl;
    
    myfile << "observed distance can = ";
    for(int i = 0; i<gen_source_can.length(); i++)
    {
      int current_gen_source = gen_source_can[i];
      if(current_gen_source != -1)
      {
        myfile << gen_matrix_can(i, current_gen_source) << " ";
      }
      else
      {
        myfile << "-1 ";
      }
    }
    myfile << std::endl;
    
    IntegerVector sample_times_cur = data["sample_times"];
    IntegerVector sample_times_can = data_can["sample_times"];
    myfile << "sample times cur = ";
    for(int i = 0; i<sample_times_cur.length(); i++)
    {
      int current_gen_source = gen_source_cur[i];
      if(current_gen_source != -1)
      {
        int time_diff = sample_times_cur[i] - sample_times_cur[current_gen_source];
        myfile << time_diff << " ";
      }
      else
      {
        myfile << "-1 ";
      }
    }
    myfile << std::endl;
    
    myfile << "sample times can = ";
    for(int i = 0; i<sample_times_can.length(); i++)
    {
      int current_gen_source = gen_source_can[i];
      if(current_gen_source != -1)
      {
        int time_diff = sample_times_can[i] - sample_times_can[current_gen_source];
        myfile << time_diff << " ";
      }
      else
      {
        myfile << "-1 ";
      }
    }
    myfile << std::endl;
  }
  
  
  if(log(U) < loglik_can - loglik_cur + log_prop_ratio)
  {
    //data = data_can;
    loglik_cur = loglik_can;
    nacc_move++;
    if(print_details) Rcout << " - accepted." << std::endl;
    if(write_details_to_file) myfile << " - accepted." << std::endl << std::endl;
    return data_can;
  }
  else
  {
    if(print_details) Rcout << " - rejected." << std::endl;
    if(write_details_to_file) myfile << " - rejected." << std::endl << std::endl;
    return data;
  }
  
}


// Update all sources without changing the time of colonisation
List UpdateAllSources(List &data, NumericVector parameters, double &loglik_cur, double w, IntegerVector Vq, IntegerVector Va, int &nacc_move, std::ofstream &myfile)
{
  bool print_details = false;
  bool write_details_to_file = true;
  
  IntegerVector t_a = data["t_a"];
  IntegerVector t_d = data["t_d"];
  IntegerVector t_c = data["t_c"];
  IntegerVector source = data["source"];
  
  List data_can = clone(data);
  //IntegerVector t_c_can = data_can["t_c"];
  IntegerVector source_can = data_can["source"];
  
  IntegerVector colonised_population = union_(Vq,Va);
  double log_prop_ratio = 0.0;
  for(int i = 0; i<colonised_population.length(); i++)
  {
    int target = colonised_population[i];
    if(source[target] != -1)
    {
      IntegerVector possible_infectors = ReturnPossibleInfectors(data, target);
      int current_source = source[target];
      IntegerVector current_source_loc = WhichVec(current_source, possible_infectors);
      if(current_source_loc.length()>0)
      {
        possible_infectors = RemoveElement(current_source_loc[0], possible_infectors);
        if(possible_infectors.length()==0)
        {
          source_can[target] = current_source;
        }
        else
        {
          source_can[target] = SampleVector(possible_infectors);
        }
      }
    }
  }
  
  


  
  bool impute_distances = true;
  
  if(data.containsElementNamed("genetic_ids") && impute_distances)
  {
    //data_can = UpdateImputedNodesMove3(data, data_can, parameters, target);
    data_can = UpdateImputedNodesMoveAll(data, data_can, parameters);
    double genetic_contribution = data_can["genetic_contribution"];
    if(genetic_contribution < -500000 || genetic_contribution > 5000000)
    {
      Rcout << "Genetic contribution = " << genetic_contribution << std::endl;
      stop("problematic genetic contribution");
    }
    log_prop_ratio += genetic_contribution;
  }
  
  
  double loglik_can = CalculateLogLikelihood(data_can, parameters);
  double U = R::runif(0.0,1.0);
  if(print_details)
  {
    Rcout << "UPDATE ALL SOURCES" << std::endl;
    Rcout << "Loglik can = " << loglik_can << ", loglik cur = " << loglik_cur << ", log prop ratio = " << log_prop_ratio;
  }
  
  if(write_details_to_file)
  {
    myfile << "UPDATE ALL SOURCES" << std::endl;
    myfile << "Loglik can = " << loglik_can << ", loglik cur = " << loglik_cur << ", log prop ratio = " << log_prop_ratio << std::endl;
    
    double transmission_ll_cur = CalculateTransmissionLikelihood(data, parameters);
    double transmission_ll_can = CalculateTransmissionLikelihood(data_can, parameters);
    
    double screening_ll_cur = CalculateScreeningLikelihood(data, parameters);
    double screening_ll_can = CalculateScreeningLikelihood(data_can, parameters);
    
    double import_ll_cur = CalculateImportationLikelihood(data, parameters);
    double import_ll_can = CalculateImportationLikelihood(data_can, parameters);
    
    double genetic_ll_cur = CalculateGeneticLikelihood(data, parameters);
    double genetic_ll_can = CalculateGeneticLikelihood(data_can, parameters);
    
    myfile << "screen ll cur = " << screening_ll_cur << ", import ll cur = " << import_ll_cur << ", transmission ll cur = "
           << transmission_ll_cur << ", genetic ll cur = " << genetic_ll_cur << std::endl;
    myfile << "screen ll can = " << screening_ll_can << ", import ll can = " << import_ll_can << ", transmission ll can = "
           << transmission_ll_can << ", genetic ll can = " << genetic_ll_can << std::endl;
    
    IntegerVector gen_source_cur = ReturnGenSourceVector(data);
    IntegerVector gen_source_can = ReturnGenSourceVector(data_can);
    
    myfile << "       Gen source cur = ";
    
    for(int i = 0; i<gen_source_cur.length(); i++)
    {
      int current_gen_source = gen_source_cur[i];
      if(current_gen_source == -1)
      {
        myfile << current_gen_source << " ";
      }
      else if(current_gen_source < 10)
      {
        myfile << " " << current_gen_source << " ";
      }
      else
      {
        myfile << current_gen_source << " ";
      }
    }
    myfile << std::endl;
    
    
    myfile << "       Gen source can = ";
    for(int i = 0; i<gen_source_can.length(); i++)
    {
      int current_gen_source = gen_source_can[i];
      if(current_gen_source == -1)
      {
        myfile << current_gen_source << " ";
      }
      else if(current_gen_source < 10)
      {
        myfile << " " << current_gen_source << " ";
      }
      else
      {
        myfile << current_gen_source << " ";
      }
    }
    myfile << std::endl;
    
    IntegerMatrix gen_matrix_cur = data["genetic_matrix"];
    IntegerMatrix gen_matrix_can = data_can["genetic_matrix"];
    
    myfile << "observed distance cur = ";
    for(int i = 0; i<gen_source_cur.length(); i++)
    {
      int current_gen_source = gen_source_cur[i];
      int distance = gen_matrix_cur(i, current_gen_source);
      if(current_gen_source != -1)
      {
        if(distance < 10)
        {
          myfile << " " << distance << " ";
        }
        else
        {
          myfile << distance << " ";
        }
      }
      else
      {
        myfile << "-1 "; 
      }
    }
    myfile << std::endl;
    
    myfile << "observed distance can = ";
    for(int i = 0; i<gen_source_can.length(); i++)
    {
      int current_gen_source = gen_source_can[i];
      int distance = gen_matrix_can(i, current_gen_source);
      if(current_gen_source != -1)
      {
        if(distance < 10)
        {
          myfile << " " << distance << " ";
        }
        else
        {
          myfile << distance << " ";
        }
      }
      else
      {
        myfile << "-1 "; 
      }
    }
    myfile << std::endl;
    
    IntegerVector sample_times_cur = data["sample_times"];
    IntegerVector sample_times_can = data_can["sample_times"];
    myfile << "     sample times cur = ";
    for(int i = 0; i<sample_times_cur.length(); i++)
    {
      int current_gen_source = gen_source_cur[i];
      if(current_gen_source != -1)
      {
        int time_diff = sample_times_cur[i] - sample_times_cur[current_gen_source];
        if(time_diff < 10)
        {
          myfile << " " << time_diff << " ";
        }
        else
        {
          myfile << time_diff << " ";
        }
        
      }
      else
      {
        myfile << "-1 ";
      }
    }
    myfile << std::endl;
    
    myfile << "     sample times can = ";
    for(int i = 0; i<sample_times_can.length(); i++)
    {
      int current_gen_source = gen_source_can[i];
      if(current_gen_source != -1)
      {
        int time_diff = sample_times_can[i] - sample_times_can[current_gen_source];
        if(time_diff < 10)
        {
          myfile << " " << time_diff << " ";
        }
        else
        {
          myfile << time_diff << " ";
        }
      }
      else
      {
        myfile << "-1 ";
      }
    }
    myfile << std::endl;
  }
  
  
  if(log(U) < loglik_can - loglik_cur + log_prop_ratio)
  {
    //data = data_can;
    loglik_cur = loglik_can;
    nacc_move++;
    if(print_details) Rcout << " - accepted." << std::endl;
    if(write_details_to_file) myfile << " - accepted." << std::endl << std::endl;
    return data_can;
  }
  else
  {
    if(print_details) Rcout << " - rejected." << std::endl;
    if(write_details_to_file) myfile << " - rejected." << std::endl << std::endl;
    return data;
  }
  
}


// Update random number of  sources without changing the time of colonisation
List UpdateRandomSources(List &data, NumericVector parameters, double &loglik_cur, double w, IntegerVector Vq, IntegerVector Va, int &nacc_move, std::ofstream &myfile)
{
  bool print_details = false;
  bool write_details_to_file = true;
  
  IntegerVector t_a = data["t_a"];
  IntegerVector t_d = data["t_d"];
  IntegerVector t_c = data["t_c"];
  IntegerVector source = data["source"];
  
  List data_can = clone(data);
  //IntegerVector t_c_can = data_can["t_c"];
  IntegerVector source_can = data_can["source"];
  
  IntegerVector colonised_population = union_(Vq,Va);
  
  
  int num_to_sample = SampleVector(seq(1,colonised_population.length()));
  IntegerVector targets_to_sample = SampleMultipleVector(colonised_population, num_to_sample);
  double log_prop_ratio = 0.0;
  for(int i = 0; i<targets_to_sample.length(); i++)
  {
    int target = colonised_population[i];
    if(source[target] != -1)
    {
      IntegerVector possible_infectors = ReturnPossibleInfectors(data, target);
      int current_source = source[target];
      IntegerVector current_source_loc = WhichVec(current_source, possible_infectors);
      if(current_source_loc.length()>0)
      {
        possible_infectors = RemoveElement(current_source_loc[0], possible_infectors);
        if(possible_infectors.length()==0)
        {
          source_can[target] = current_source;
        }
        else
        {
          source_can[target] = SampleVector(possible_infectors);
        }
      }
    }
  }
  
  
  
  
  
  bool impute_distances = true;
  
  if(data.containsElementNamed("genetic_ids") && impute_distances)
  {
    //data_can = UpdateImputedNodesMove3(data, data_can, parameters, target);
    data_can = UpdateImputedNodesMoveAll(data, data_can, parameters);
    double genetic_contribution = data_can["genetic_contribution"];
    if(genetic_contribution < -500000 || genetic_contribution > 5000000)
    {
      Rcout << "Genetic contribution = " << genetic_contribution << std::endl;
      stop("problematic genetic contribution");
    }
    log_prop_ratio += genetic_contribution;
  }
  
  
  double loglik_can = CalculateLogLikelihood(data_can, parameters);
  double U = R::runif(0.0,1.0);
  if(print_details)
  {
    Rcout << "UPDATE "<< num_to_sample << " SOURCES" << std::endl;
    Rcout << "SOURCES = " << targets_to_sample << std::endl;
    Rcout << "Loglik can = " << loglik_can << ", loglik cur = " << loglik_cur << ", log prop ratio = " << log_prop_ratio;
  }
  
  if(write_details_to_file)
  {
    myfile << "UPDATE "<< num_to_sample << " SOURCES" << std::endl;
    myfile << "SOURCES = " << targets_to_sample << std::endl;
    myfile << "Loglik can = " << loglik_can << ", loglik cur = " << loglik_cur << ", log prop ratio = " << log_prop_ratio << std::endl;
    
    double transmission_ll_cur = CalculateTransmissionLikelihood(data, parameters);
    double transmission_ll_can = CalculateTransmissionLikelihood(data_can, parameters);
    
    double screening_ll_cur = CalculateScreeningLikelihood(data, parameters);
    double screening_ll_can = CalculateScreeningLikelihood(data_can, parameters);
    
    double import_ll_cur = CalculateImportationLikelihood(data, parameters);
    double import_ll_can = CalculateImportationLikelihood(data_can, parameters);
    
    double genetic_ll_cur = CalculateGeneticLikelihood(data, parameters);
    double genetic_ll_can = CalculateGeneticLikelihood(data_can, parameters);
    
    myfile << "screen ll cur = " << screening_ll_cur << ", import ll cur = " << import_ll_cur << ", transmission ll cur = "
           << transmission_ll_cur << ", genetic ll cur = " << genetic_ll_cur << std::endl;
    myfile << "screen ll can = " << screening_ll_can << ", import ll can = " << import_ll_can << ", transmission ll can = "
           << transmission_ll_can << ", genetic ll can = " << genetic_ll_can << std::endl;
    
    IntegerVector gen_source_cur = ReturnGenSourceVector(data);
    IntegerVector gen_source_can = ReturnGenSourceVector(data_can);
    
    myfile << "       Gen source cur = ";
    
    for(int i = 0; i<gen_source_cur.length(); i++)
    {
      int current_gen_source = gen_source_cur[i];
      if(current_gen_source == -1)
      {
        myfile << current_gen_source << " ";
      }
      else if(current_gen_source < 10)
      {
        myfile << " " << current_gen_source << " ";
      }
      else
      {
        myfile << current_gen_source << " ";
      }
    }
    myfile << std::endl;
    
    
    myfile << "       Gen source can = ";
    for(int i = 0; i<gen_source_can.length(); i++)
    {
      int current_gen_source = gen_source_can[i];
      if(current_gen_source == -1)
      {
        myfile << current_gen_source << " ";
      }
      else if(current_gen_source < 10)
      {
        myfile << " " << current_gen_source << " ";
      }
      else
      {
        myfile << current_gen_source << " ";
      }
    }
    myfile << std::endl;
    
    IntegerMatrix gen_matrix_cur = data["genetic_matrix"];
    IntegerMatrix gen_matrix_can = data_can["genetic_matrix"];
    
    myfile << "observed distance cur = ";
    for(int i = 0; i<gen_source_cur.length(); i++)
    {
      int current_gen_source = gen_source_cur[i];
      int distance = gen_matrix_cur(i, current_gen_source);
      if(current_gen_source != -1)
      {
        if(distance < 10)
        {
          myfile << " " << distance << " ";
        }
        else
        {
          myfile << distance << " ";
        }
      }
      else
      {
        myfile << "-1 "; 
      }
    }
    myfile << std::endl;
    
    myfile << "observed distance can = ";
    for(int i = 0; i<gen_source_can.length(); i++)
    {
      int current_gen_source = gen_source_can[i];
      int distance = gen_matrix_can(i, current_gen_source);
      if(current_gen_source != -1)
      {
        if(distance < 10)
        {
          myfile << " " << distance << " ";
        }
        else
        {
          myfile << distance << " ";
        }
      }
      else
      {
        myfile << "-1 "; 
      }
    }
    myfile << std::endl;
    
    IntegerVector sample_times_cur = data["sample_times"];
    IntegerVector sample_times_can = data_can["sample_times"];
    myfile << "     sample times cur = ";
    for(int i = 0; i<sample_times_cur.length(); i++)
    {
      int current_gen_source = gen_source_cur[i];
      if(current_gen_source != -1)
      {
        int time_diff = sample_times_cur[i] - sample_times_cur[current_gen_source];
        if(time_diff < 10)
        {
          myfile << " " << time_diff << " ";
        }
        else
        {
          myfile << time_diff << " ";
        }
        
      }
      else
      {
        myfile << "-1 ";
      }
    }
    myfile << std::endl;
    
    myfile << "     sample times can = ";
    for(int i = 0; i<sample_times_can.length(); i++)
    {
      int current_gen_source = gen_source_can[i];
      if(current_gen_source != -1)
      {
        int time_diff = sample_times_can[i] - sample_times_can[current_gen_source];
        if(time_diff < 10)
        {
          myfile << " " << time_diff << " ";
        }
        else
        {
          myfile << time_diff << " ";
        }
      }
      else
      {
        myfile << "-1 ";
      }
    }
    myfile << std::endl;
  }
  
  
  if(log(U) < loglik_can - loglik_cur + log_prop_ratio)
  {
    //data = data_can;
    loglik_cur = loglik_can;
    nacc_move++;
    if(print_details) Rcout << " - accepted." << std::endl;
    if(write_details_to_file) myfile << " - accepted." << std::endl << std::endl;
    return data_can;
  }
  else
  {
    if(print_details) Rcout << " - rejected." << std::endl;
    if(write_details_to_file) myfile << " - rejected." << std::endl << std::endl;
    return data;
  }
  
}


// Swap a target with their source
List SwapTargetOffspring(List &data, NumericVector parameters, double &loglik_cur, double w, IntegerVector Vq, IntegerVector Va, int &nacc_move, std::ofstream &myfile)
{
  bool print_details = false;
  bool write_details_to_file = true;
  
  IntegerVector t_a = data["t_a"];
  IntegerVector t_d = data["t_d"];
  IntegerVector t_c = data["t_c"];
  IntegerVector source = data["source"];
  
  List data_can = clone(data);
  IntegerVector t_c_can = data_can["t_c"];
  IntegerVector source_can = data_can["source"];
  
  IntegerVector colonised_population = union_(Vq,Va);
  int target = SampleVector(colonised_population);
  int num_patients = data["num_patients"];

  //int last_day = ReturnLastDay(data, target);
  
  bool print_debug = false;
  //if(target == 59 || target == 62 || target == 56) print_debug = true;
  
  double log_prop_ratio = 0.0;
  
  int target_source = source[target];
  
  if(print_debug)
  {
    Rcout << "Target = " << target << ", source = " << target_source << " - attempt to swap" << std::endl;
  }
  
  // remove if the source is an importation or has a source of hcw
  if(target_source == -1) return data;
  if(target_source >= num_patients) return data;
  
  t_c_can[target] = t_c[target_source];
  t_c_can[target_source] = t_c[target];
  

  
  source_can[target] = source[target_source];
  source_can[target_source] = target;
  
  
  // Now check if the colonisation time of the source is before the admission of the target
  if(t_c_can[target] < t_a[target]) return data;
  
  // Check if the colonisation time of the source is after the birth of another offspring or a swab
  IntegerVector source_infections = WhichVec(target_source, source);
  for(int i = 0; i<source_infections.length(); i++)
  {
    int current_infection_target = source_infections[i];
    if(current_infection_target != target)
    {
      int current_infection_target_col_time = t_c[current_infection_target];
      if(current_infection_target_col_time <= t_c_can[target_source])
      {
        return data;
      }
    }
  }
  
  IntegerMatrix screening_matrix = data["screening_matrix"];
  for(int i = 0; i<screening_matrix.ncol(); i++)
  {
    if(screening_matrix(target_source,i)==1)
    {
      // Source of the target has a positive swab, check if before the new col time
      if(t_c_can[target_source] > i) return data;
    }
  }
  
  
  
  bool impute_distances = true;
  
  if(data.containsElementNamed("genetic_ids") && impute_distances)
  {
    //data_can = UpdateImputedNodesMove3(data, data_can, parameters, target);
    data_can = UpdateImputedNodesMoveAll(data, data_can, parameters);
    double genetic_contribution = data_can["genetic_contribution"];
    if(genetic_contribution < -500000 || genetic_contribution > 5000000)
    {
      Rcout << "Genetic contribution = " << genetic_contribution << std::endl;
      stop("problematic genetic contribution");
    }
    log_prop_ratio += genetic_contribution;
  }
  
  
  double loglik_can = CalculateLogLikelihood(data_can, parameters);
  double U = R::runif(0.0,1.0);
  if(print_details)
  {
    Rcout << "SWAP SOURCE OFFSPRING ID = " << target << " from source " << source[target] << " to source " << source_can[target] << std::endl;
    Rcout << "Loglik can = " << loglik_can << ", loglik cur = " << loglik_cur << ", log prop ratio = " << log_prop_ratio;
  }
  
  if(write_details_to_file)
  {
    myfile << "SWAP SOURCE OFFSPRING ID = " << target << " from source " << source[target] << " to source " << source_can[target] << std::endl;
    myfile << "Loglik can = " << loglik_can << ", loglik cur = " << loglik_cur << ", log prop ratio = " << log_prop_ratio << std::endl;
    
    double transmission_ll_cur = CalculateTransmissionLikelihood(data, parameters);
    double transmission_ll_can = CalculateTransmissionLikelihood(data_can, parameters);
    
    double screening_ll_cur = CalculateScreeningLikelihood(data, parameters);
    double screening_ll_can = CalculateScreeningLikelihood(data_can, parameters);
    
    double import_ll_cur = CalculateImportationLikelihood(data, parameters);
    double import_ll_can = CalculateImportationLikelihood(data_can, parameters);
    
    double genetic_ll_cur = CalculateGeneticLikelihood(data, parameters);
    double genetic_ll_can = CalculateGeneticLikelihood(data_can, parameters);
    
    myfile << "screen ll cur = " << screening_ll_cur << ", import ll cur = " << import_ll_cur << ", transmission ll cur = "
           << transmission_ll_cur << ", genetic ll cur = " << genetic_ll_cur << std::endl;
    myfile << "screen ll can = " << screening_ll_can << ", import ll can = " << import_ll_can << ", transmission ll can = "
           << transmission_ll_can << ", genetic ll can = " << genetic_ll_can << std::endl;
    
    IntegerVector gen_source_cur = ReturnGenSourceVector(data);
    IntegerVector gen_source_can = ReturnGenSourceVector(data_can);
    
    myfile << "Gen source cur = " << gen_source_cur << std::endl;
    myfile << "Gen source can = " << gen_source_can << std::endl;
    
    IntegerMatrix gen_matrix_cur = data["genetic_matrix"];
    IntegerMatrix gen_matrix_can = data_can["genetic_matrix"];
    
    myfile << "observed distance cur = ";
    for(int i = 0; i<gen_source_cur.length(); i++)
    {
      int current_gen_source = gen_source_cur[i];
      if(current_gen_source != -1)
      {
        myfile << gen_matrix_cur(i, current_gen_source) << " ";
      }
      else
      {
        myfile << "-1 "; 
      }
      
    }
    myfile << std::endl;
    
    myfile << "observed distance can = ";
    for(int i = 0; i<gen_source_can.length(); i++)
    {
      int current_gen_source = gen_source_can[i];
      if(current_gen_source != -1)
      {
        myfile << gen_matrix_can(i, current_gen_source) << " ";
      }
      else
      {
        myfile << "-1 ";
      }
    }
    myfile << std::endl;
    
    IntegerVector sample_times_cur = data["sample_times"];
    IntegerVector sample_times_can = data_can["sample_times"];
    myfile << "sample times cur = ";
    for(int i = 0; i<sample_times_cur.length(); i++)
    {
      int current_gen_source = gen_source_cur[i];
      if(current_gen_source != -1)
      {
        int time_diff = sample_times_cur[i] - sample_times_cur[current_gen_source];
        myfile << time_diff << " ";
      }
      else
      {
        myfile << "-1 ";
      }
    }
    myfile << std::endl;
    
    myfile << "sample times can = ";
    for(int i = 0; i<sample_times_can.length(); i++)
    {
      int current_gen_source = gen_source_can[i];
      if(current_gen_source != -1)
      {
        int time_diff = sample_times_can[i] - sample_times_can[current_gen_source];
        myfile << time_diff << " ";
      }
      else
      {
        myfile << "-1 ";
      }
    }
    myfile << std::endl;
  }
  
  
  if(log(U) < loglik_can - loglik_cur + log_prop_ratio)
  {
    //data = data_can;
    loglik_cur = loglik_can;
    nacc_move++;
    if(print_details) Rcout << " - accepted." << std::endl;
    if(write_details_to_file) myfile << " - accepted." << std::endl << std::endl;
    return data_can;
  }
  else
  {
    if(print_details) Rcout << " - rejected." << std::endl;
    if(write_details_to_file) myfile << " - rejected." << std::endl << std::endl;
    return data;
  }
  
}

void MoveColonisationTime(List &data, NumericVector parameters, double &loglik_cur, double w, IntegerVector Vq, IntegerVector Va, int &nacc_move)
{
  bool print_details = false;
  
  IntegerVector t_a = data["t_a"];
  IntegerVector t_d = data["t_d"];
  IntegerVector t_c = data["t_c"];
  IntegerVector source = data["source"];
  
  List data_can = clone(data);
  IntegerVector t_c_can = data_can["t_c"];
  IntegerVector source_can = data_can["source"];
  
  IntegerVector colonised_population = union_(Vq,Va);
  int target = SampleVector(colonised_population);
  //Rcout << "Colonised population = " << colonised_population << std::endl;
  int last_day = ReturnLastDay(data, target);
  //Rcout << "t_a = " << t_a[target] << "Last day = " << last_day << std::endl;
  
  
  double log_prop_ratio = 0.0;
  double U = R::runif(0.0,1.0);
  if(U < w)
  {
    // Propose the target is an importation
    t_c_can[target] = t_a[target];
    source_can[target] = -1;
    if(source[target]==-1)
    {
      // Importation -> Importation
      // Do nothing, proposal ratio = 1
    }
    else
    {
      // Acquisition -> Importation
      log_prop_ratio = log(1-w) - log(w) - log(last_day-t_a[target]+1) - log(total_col_pop(data, t_c[target]));
    }
  }
  else
  {
    // Propose the target is an acquisition
    if(last_day < t_a[target]) return; // The target cannot be an acquisition
    t_c_can[target] = SampleVector(seq(t_a[target],last_day));
    IntegerVector possible_infectors = ReturnPossibleInfectors(data_can, target);
    //Rcout << "Possible infectors " << possible_infectors << std::endl;
    if(possible_infectors.length()==0) return;
    source_can[target] = SampleVector(possible_infectors);
    if(source[target] != -1)
    {
      // Acquisition -> Acquisition
      log_prop_ratio = log(total_col_pop(data_can, t_c_can[target])) - log(total_col_pop(data, t_c[target]));
    }
    else
    {
      // Importation -> Acquisition
      log_prop_ratio = log(w) + log(last_day - t_a[target] + 1) - log(1-w) + log(total_col_pop(data_can, t_c_can[target]));
    }
  }
  
  //if(t_c[target] == t_c_can[target] && source[target] == source_can[target]) return;
  
  
  double loglik_can = CalculateLogLikelihood(data_can, parameters);
  U = R::runif(0.0,1.0);
  if(print_details)
  {
    Rcout << "  MOVE ID = " << target << " from time t = " << t_c[target] << " to t = " << t_c_can[target]
          << " from source " << source[target] << " to source " << source_can[target] << std::endl;
    Rcout << "Loglik can = " << loglik_can << ", loglik cur = " << loglik_cur << ", log prop ratio = " << log_prop_ratio;
    Rcout << std::endl;
  }
  if(log(U) < loglik_can - loglik_cur + log_prop_ratio)
  {
    data = data_can;
    loglik_cur = loglik_can;
    nacc_move++;
    if(print_details) Rcout << " - accepted." << std::endl;
  }
  else
  {
    if(print_details) Rcout << " - rejected." << std::endl;
  }
}

void MoveColonisationTime2(List &data, NumericVector parameters, double &loglik_cur, double w, IntegerVector Vq, IntegerVector Va, int &nacc_move)
{
  bool print_details = false;
  
  IntegerVector t_a = data["t_a"];
  IntegerVector t_d = data["t_d"];
  IntegerVector t_c = data["t_c"];
  IntegerVector source = data["source"];
  
  List data_can = clone(data);
  IntegerVector t_c_can = data_can["t_c"];
  IntegerVector source_can = data_can["source"];
  
  IntegerVector colonised_population = union_(Vq,Va);
  int target = SampleVector(colonised_population);
  int num_patients = data["num_patients"];
  //Rcout << "Colonised population = " << colonised_population << std::endl;
  int last_day = ReturnLastDay(data, target);
  //Rcout << "t_a = " << t_a[target] << "Last day = " << last_day << std::endl;
  
  
  double log_prop_ratio = 0.0;
  double U = R::runif(0.0,1.0);
  double v = 0.2;
  if(U < w)
  {
    // Propose the target is an importation
    t_c_can[target] = t_a[target];
    source_can[target] = -1;
    if(source[target]==-1)
    {
      // Importation -> Importation
      // Do nothing, proposal ratio = 1
    }
    else if(source[target] < num_patients)
    {
      // Patient Acquisition -> Importation
      log_prop_ratio = log(1-w) + log(1-v) - log(w) - log(last_day-t_a[target]+1) - log(col_pop(data, t_c[target]));
    }
    else
    {
      // HCW Acquisition -> Importation
      log_prop_ratio = log(1-w) + log(v) - log(w) - log(last_day-t_a[target]+1) - log(hcw_col_pop(data, t_c[target]));
    }
  }
  else
  {
    // Propose the target is an acquisition
    if(last_day < t_a[target]) return; // The target cannot be an acquisition
    t_c_can[target] = SampleVector(seq(t_a[target],last_day));
    
    U = R::runif(0.0,1.0);
    if(U < v)
    {
      // We are proposing an acquisition from a healthcare worker
      IntegerVector possible_infectors = ReturnPossibleInfectors_hcw(data_can, target);
      if(possible_infectors.length()==0) return;
      source_can[target] = SampleVector(possible_infectors);
      if(source[target] == -1)
      {
        // Importation -> HCW Acquisition
        log_prop_ratio = log(w) + log(last_day-t_a[target]+1) + log(hcw_col_pop(data_can,t_c_can[target])) - log(1-w) - log(v);
      }
      else if(source[target] < num_patients)
      {
        // Patient Acquisition -> HCW Acquisition
        log_prop_ratio = log(1-v) + log(hcw_col_pop(data_can,t_c_can[target]) - log(v) - log(col_pop(data,t_c[target])));
      }
      else
      {
        // HCW Acquisition -> HCW Acquisition
        log_prop_ratio = log(hcw_col_pop(data_can,t_c_can[target])) - log(hcw_col_pop(data,t_c[target]));			
      }
      if(hcw_col_pop(data_can,t_c_can[target]) != possible_infectors.length())
      {
        stop("inconsistency between possible infectors and the hcw_col_pop");
      }
    }
    else
    {
      // We are proposing an acquisition from a patient
      IntegerVector possible_infectors = ReturnPossibleInfectors_patient(data_can, target);
      if(possible_infectors.length()==0) return;
      source_can[target] = SampleVector(possible_infectors);
      if(source[target] == -1)
      {
        // Importation -> Patient Acquisition
        log_prop_ratio = log(w) + log(last_day-t_a[target]+1) + log(col_pop(data_can,t_c_can[target])) - log(1-w) - log(1-v);
      }
      else if(source[target] < num_patients)
      {
        // Patient Acquisition -> Patient Acquisition
        log_prop_ratio = log(col_pop(data_can,t_c_can[target])) - log(col_pop(data,t_c[target]));	
      }
      else
      {
        // HCW Acquisition -> Patient Acquisition
        log_prop_ratio = log(v) + log(col_pop(data_can,t_c_can[target])) - log(1-v) - log(hcw_col_pop(data,t_c[target]));
      }
      if(col_pop(data_can,t_c_can[target]) != possible_infectors.length())
      {
        stop("inconsistency between possible infectors and the col_pop");
      }
    }
    
  }
  
  
  double loglik_can = CalculateLogLikelihood(data_can, parameters);
  U = R::runif(0.0,1.0);
  if(print_details)
  {
    Rcout << std::endl;
    Rcout << "  MOVE ID = " << target << " from time t = " << t_c[target] << " to t = " << t_c_can[target]
          << " from source " << source[target] << " to source " << source_can[target] << std::endl;
    Rcout << "Loglik can = " << loglik_can << ", loglik cur = " << loglik_cur << ", log prop ratio = " << log_prop_ratio;
  }
  if(log(U) < loglik_can - loglik_cur + log_prop_ratio)
  {
    data = data_can;
    loglik_cur = loglik_can;
    nacc_move++;
    if(print_details) Rcout << " - accepted." << std::endl;
  }
  else
  {
    if(print_details) Rcout << " - rejected." << std::endl;
  }
}


void AddColonisationTime(List &data, NumericVector parameters, double &loglik_cur, double w, IntegerVector Vs, IntegerVector &Va, int &nacc_add)
{
  bool print_details = false;
  IntegerVector t_a = data["t_a"];
  IntegerVector t_c = data["t_c"];
  IntegerVector t_d = data["t_d"];
  IntegerVector source = data["source"];
  
  List data_can = clone(data);
  IntegerVector t_c_can = data_can["t_c"];
  IntegerVector source_can = data_can["source"];
  
  IntegerVector susceptible_population = setdiff(Vs,Va);
  if(susceptible_population.length()==0) return;
  IntegerVector nadd0 = ReturnNadd0(data, Va);
  int target = SampleVector(susceptible_population);
  
  double log_prop_ratio = 0;
  
  double U = R::runif(0.0,1.0);
  if(U < w)
  {
    // Propose an importation
    t_c_can[target] = t_a[target];
    source_can[target] = -1;
    
    log_prop_ratio = log(susceptible_population.length()) - log(w) - log(1+nadd0.length());
    //log_prop_ratio = log(susceptible_population.length()) - log(w) - log(1+Va.length());
  }
  else
  {
    // Propose an acquisition
    t_c_can[target] = SampleVector(seq(t_a[target],t_d[target]));
    IntegerVector possible_infectors = ReturnPossibleInfectors(data_can, target);
    if(possible_infectors.length()==0) return;
    source_can[target] = SampleVector(possible_infectors);
    
    log_prop_ratio = log(total_col_pop(data_can, t_c_can[target])) + log(susceptible_population.length()) + 
      log(t_d[target] - t_a[target] + 1) - log(1-w) - log(1+nadd0.length());
  }
  
  //if(!CheckData(data_can)) stop("error with data with the adding step!");
  double loglik_can = CalculateLogLikelihood(data_can, parameters);
  U = R::runif(0.0,1.0);
  
  if(print_details) Rcout << "ADDING ID " << target <<  ": Loglik can = " << loglik_can << ", loglik cur = " << loglik_cur << ", log prop ratio = " << log_prop_ratio;
  if(log(U) < loglik_can - loglik_cur + log_prop_ratio)
  {
    data = data_can;
    loglik_cur = loglik_can;
    Va.push_back(target);
    nacc_add++;
    if(print_details) Rcout << " - accepted. Number added = " << Va.length() << std::endl;
  }
  else
  {
    if(print_details) Rcout << " - rejected." << std::endl;
  }
  
}

void RemoveColonisationTime(List &data, NumericVector parameters, double &loglik_cur, double w, IntegerVector Vs, IntegerVector &Va, int &nacc_remove)
{
  bool print_details = false;
  IntegerVector t_a = data["t_a"];
  IntegerVector t_c = data["t_c"];
  IntegerVector t_d = data["t_d"];
  IntegerVector source = data["source"];
  
  List data_can = clone(data);
  IntegerVector t_c_can = data_can["t_c"];
  IntegerVector source_can = data_can["source"];
  
  IntegerVector susceptible_population = setdiff(Vs,Va);
  IntegerVector nadd0 = ReturnNadd0(data, Va);
  if(nadd0.length()==0) return;
  int target = SampleVector(nadd0);
  double log_prop_ratio = 0;
  
  source_can[target] = -2;
  t_c_can[target] = -1;
  
  if(source[target]==-1)
  {
    // Removing an importation
    log_prop_ratio = log(nadd0.length()) + log(w) - log(susceptible_population.length()+1);
    //log_prop_ratio = log(Va.length()) + log(w) - log(susceptible_population.length()+1);
  }
  else
  {
    // Removing an acquisition
    log_prop_ratio = log(nadd0.length()) + log(1-w) - log(t_d[target]-t_a[target]+1) - log(susceptible_population.length()+1) - log(total_col_pop(data,t_c[target]));
    //log_prop_ratio = log(Va.length()) + log(1-w) - log(t_d[target]-t_a[target]+1) - log(susceptible_population.length()+1) - log(col_pop(data,t_c[target]));
    
  }
  
  //if(!CheckData(data_can)) stop("error with data in the remove step!");
  //Rcout << "Before likelihood" << std::endl;
  double loglik_can = CalculateLogLikelihood(data_can, parameters);
  //Rcout << "After likelihood" << std::endl;
  double U = R::runif(0.0,1.0);
  if(print_details) Rcout << "Removing ID " << target <<  ": Loglik can = " << loglik_can << ", loglik cur = " << loglik_cur << ", log prop ratio = " << log_prop_ratio;
  if(log(U) < loglik_can - loglik_cur + log_prop_ratio)
  {
    data = data_can;
    loglik_cur = loglik_can;
    int Va_loc = WhichVec(target, Va)[0];
    Va = RemoveElement(Va_loc, Va);
    if(print_details) Rcout << " - accepted. Number added = " << Va.length() << std::endl;
    nacc_remove++;
  }
  else
  {
    if(print_details) Rcout << " - rejected." << std::endl;
  }
}


bool DoesTargetHaveSequenceGreaterThanTime(List data, int time, int target, int current_variant)
{
  IntegerVector genetic_ids = data["genetic_ids"];
  IntegerVector sample_times = data["sample_times"];
  IntegerVector variant_numbers = data["variant_numbers"];
  IntegerVector variant_idx = WhichVec(current_variant, variant_numbers);
  IntegerVector target_genetic_idx = WhichVec(target, genetic_ids);
  IntegerVector target_variant_genetic_idx = intersect(target_genetic_idx, variant_idx);
  
  IntegerVector target_sample_times = as<IntegerVector>(sample_times[target_variant_genetic_idx]);
  for(int i = 0; i<target_sample_times.length(); i++)
  {
    int current_time = target_sample_times[i];
    if(current_time > time)
    {
      return true;
    }
  }
  return false;
}



// [[Rcpp::export]]
void MCMC_EPI_SOURCE(List MCMC_options, IntegerVector t_a, IntegerVector t_c, IntegerVector t_d, IntegerVector source, 
                     IntegerVector hcw_ind, IntegerMatrix screening_matrix)
{
  // Initialise data
  List data = InitialiseData(t_a, t_c, t_d, source, hcw_ind, screening_matrix);
  int num_patients = data["num_patients"];
  
  // Load MCMC options
  const int iterations = MCMC_options["iterations"];
  const int num_updates = MCMC_options["num_updates"];
  NumericVector prior_parameters = MCMC_options["prior_parameters"];
  NumericVector initial_chain_state = MCMC_options["initial_chain_state"];
  NumericVector parameters = clone(initial_chain_state);
  NumericVector proposal_variance = MCMC_options["proposal_variance"];
  IntegerVector debug_flags = MCMC_options["debug_flags"];
  
  std::string output_file = MCMC_options["output_file"];
  std::string source_file = MCMC_options["source_file"];
  std::string coltime_file = MCMC_options["coltime_file"];
  
  // Check if there are healthcare workers
  bool include_hcw = true;
  if(t_c.length() == num_patients)
  {
    include_hcw = false;
    parameters[3] = 0.0;
    Rcout << "No healthcare workers present - use standard model" << std::endl;
  }
  else
  {
    Rcout << "Healthcare workers present, include beta_h term" << std::endl;
  }
  
  // Calculate other variables
  int true_positives = CalculateTruePositives(data);
  int nacc_beta_p = 0;
  int nacc_beta_h = 0;
  int nacc_add = 0;
  int nacc_move = 0;
  int nacc_remove = 0;
  double loglik = 0.0;
  IntegerVector Va;
  IntegerVector Vs = CalculateVs(data);
  IntegerVector Vq = CalculateVq(data);
  IntegerVector augmented_moves_proposed(5);
  // Write to files
  // Output file
  remove(output_file.c_str());
  std::ofstream myfile; // Define output stream
  myfile.open(output_file.c_str()); // Open file
  assert(myfile.is_open());
  PrintNumVectorToFile(myfile, parameters);
  
  // Source file
  remove(source_file.c_str());
  std::ofstream myfile2;
  myfile2.open(source_file.c_str());
  assert(myfile2.is_open());
  PrintIntVectorToFile(myfile2, data["source"]);
  
  // Col time file
  remove(coltime_file.c_str());
  std::ofstream myfile3;
  myfile3.open(coltime_file.c_str());
  assert(myfile3.is_open());
  PrintColonisationTimeSumToFile(myfile3, data);
  //PrintIntVectorToFile(myfile3, data["t_c"]);
  
  
  
  
  // Begin MCMC
  Rcout << "Begin MCMC" << std::endl;
  for(int i = 1; i < iterations; i++)
  {
    // Update z by gibbs
    if(debug_flags[0]==0)
    {
      parameters[0] = R::rbeta(prior_parameters[0] + true_positives, prior_parameters[1] + CalculateFalseNegatives(data));
    }
    
    // Update p by gibbs
    if(debug_flags[1]==0)
    {
      int import_sum = CalculateImportationSum(data);
      parameters[1] = R::rbeta(prior_parameters[2] + import_sum, prior_parameters[3] + num_patients - import_sum);
    }
    
    
    loglik = CalculateLogLikelihood(data, parameters);
    
    // Update beta_P by metropolis hastings
    if(debug_flags[2]==0)
    {
      UpdateTransmissionRate_P(data, parameters, loglik, proposal_variance[0], prior_parameters[4], nacc_beta_p);
    }
    
    // update beta_H by metropolis hastings
    if(debug_flags[3]==0 && include_hcw)
    {
      UpdateTransmissionRate_H(data, parameters, loglik, proposal_variance[1], prior_parameters[5], nacc_beta_h);
    }
    
    
    // Augmented data updates
    if(debug_flags[4]==0)
    {
      for(int update_counter = 0; update_counter < num_updates; update_counter++)
      {
        double w = 0.3;
        //double w = 0.00001;
        double move = R::runif(0.0,3.0);
        augmented_moves_proposed[floor(move)]++;
        if(floor(move) < 1)
        {
          // Move a colonisation time
          MoveColonisationTime(data, parameters, loglik, w, Vq, Va, nacc_move);
        }
        else if(floor(move) < 2)
        {
          // Add a colonisation time
          AddColonisationTime(data, parameters, loglik, w, Vs, Va, nacc_add);
        }
        else if(floor(move) < 3)
        {
          // Remove a colonisation time
          RemoveColonisationTime(data, parameters, loglik, w, Vs, Va, nacc_remove);
        }
        
        
        
        
      }
      // Write colonisation time sum to file
      PrintIntVectorToFile(myfile2, data["source"]);
      PrintColonisationTimeSumToFile(myfile3, data);
      //PrintIntVectorToFile(myfile3, data["t_c"]);
      
    }
    
    
    
    
    // Write parameters to file
    PrintNumVectorToFile(myfile, parameters);
    
    
    if(i%100==0)
    {
      Rcout << "Iteraion " << i << " completed, number added by algorithm = " << Va.length() << std::endl;
    }
    
  }
  
  // Close files
  myfile.close();
  myfile2.close();
  myfile3.close();
  
  double beta_h_prob = (double)nacc_beta_h/(double)iterations;
  double beta_p_prob = (double)nacc_beta_p/(double)iterations;
  double move_prob = (double)nacc_move/(double)(augmented_moves_proposed[0]);
  double add_prob = (double)nacc_add/(double)(augmented_moves_proposed[1]);
  double remove_prob = (double)nacc_remove/(double)(augmented_moves_proposed[2]);
  Rcout << "beta_P acceptance = " << beta_p_prob << ", beta_H acceptance = " << beta_h_prob << std::endl;
  Rcout << "Move prob = " << move_prob << ", add prob = " << add_prob << ", remove prob = " << remove_prob << std::endl;
  //Rcout << "Move prob = " << move_prob << ", add prob = " << add_prob << ", remove prob = " << remove_prob << std::endl;
}



// No source functions



double CalculateTransmissionLikelihood_NS(List data, NumericVector parameters)
{
  double loglik = 0.0;
  bool print_details = false;
  int num_patients = data["num_patients"];
  IntegerVector t_a = data["t_a"];
  IntegerVector t_c = data["t_c"];
  IntegerVector t_d = data["t_d"];
  IntegerVector source = data["source"];
  double beta_p = parameters[2];
  double beta_h = 0.0;
  if(num_patients != t_c.length())
  {
    beta_h = parameters[3];
  }
  
  for(int i = 0; i<num_patients; i++)
  {
    if(t_c[i] == -1)
    {
      // The individual avoided infection the whole time
      for(int t = t_a[i]; t<=t_d[i]; t++)
      {
        loglik -= (beta_p*col_pop(data,t) + beta_h*hcw_col_pop(data,t));
      }
    }
    else if(source[i] == -1)
    {
      // The individual is an importation, and so no contribution to the likelihood
    }
    else if(t_c[i] == t_a[i])
    {
      // The individual was infected on their first day
      loglik += log(1-exp(-(beta_p*col_pop(data,t_a[i]) + beta_h*hcw_col_pop(data,t_a[i]))));
    }
    else
    {
      // The individual was infected after the first day
      for(int t = t_a[i]; t<t_c[i]; t++)
      {
        loglik -= (beta_p*col_pop(data,t) + beta_h*hcw_col_pop(data,t));
      }
      loglik += log(1-exp(-(beta_p*col_pop(data,t_c[i]) + beta_h*hcw_col_pop(data,t_c[i]))));
    }
    if(print_details)
    {
      Rcout << "i = " << i << ", loglik = " << loglik << ", col pop = " << col_pop(data,t_c[i]) 
            << ", hcw col pop " << hcw_col_pop(data,t_c[i]) << ", source = " << source[i] << std::endl;
    }
  }
  return loglik;
}

// [[Rcpp::export]]
double CalculateLogLikelihood_NS(List data, NumericVector parameters)
{
  bool print_details = false;
  double screening_likelihood = CalculateScreeningLikelihood(data, parameters);
  double importation_likelihood = CalculateImportationLikelihood(data, parameters);
  double transmission_likelihood = CalculateTransmissionLikelihood_NS(data, parameters);
  double loglik = screening_likelihood + importation_likelihood + transmission_likelihood;
  if(print_details)
  {
    Rcout << "Screen loglik = " << screening_likelihood << ", import loglik = " << importation_likelihood << ", transmission loglik = " 
          << transmission_likelihood << std::endl;
  }
  return loglik;
}

void UpdateTransmissionRate_P_NS(List data, NumericVector &parameters, double &loglik_cur, double proposal_variance, double prior, int &nacc)
{
  bool print_details = false;
  NumericVector parameters_can = clone(parameters);
  parameters_can[2] = R::rnorm(parameters[2], proposal_variance);
  if(parameters_can[2] > 0)
  {
    double logpi_cur = loglik_cur + R::dexp(parameters[2], 1/prior, 1);
    double loglik_can = CalculateLogLikelihood_NS(data, parameters_can);
    double logpi_can = loglik_can + R::dexp(parameters_can[2], 1/prior, 1);
    double U = R::runif(0.0,1.0);
    if(print_details)
    {
      Rcout << "Beta cur = " << parameters[2] << ", beta can = " << parameters_can[2]
            << "loglik cur = " << loglik_cur << ", loglik can = " << loglik_can;
    }
    if(log(U) < logpi_can - logpi_cur)
    {
      if(print_details)
      {
        Rcout << " - accepted." << std::endl;
      }
      parameters = parameters_can;
      loglik_cur = loglik_can;
      nacc++;
    }
    else
    {
      if(print_details)
      {
        Rcout << " - rejected." << std::endl;
      }
    }
  }
}

void UpdateTransmissionRate_H_NS(List data, NumericVector &parameters, double &loglik_cur, double proposal_variance, double prior, int &nacc)
{
  bool print_details = false;
  NumericVector parameters_can = clone(parameters);
  parameters_can[3] = R::rnorm(parameters[3], proposal_variance);
  if(parameters_can[3] > 0)
  {
    double logpi_cur = loglik_cur + R::dexp(parameters[3], 1/prior, 1);
    double loglik_can = CalculateLogLikelihood_NS(data, parameters_can);
    double logpi_can = loglik_can + R::dexp(parameters_can[3], 1/prior, 1);
    double U = R::runif(0.0,1.0);
    if(print_details)
    {
      Rcout << "Beta cur = " << parameters[3] << ", beta can = " << parameters_can[3]
            << "loglik cur = " << loglik_cur << ", loglik can = " << loglik_can;
    }
    if(log(U) < logpi_can - logpi_cur)
    {
      if(print_details)
      {
        Rcout << " - accepted." << std::endl;
      }
      parameters = parameters_can;
      loglik_cur = loglik_can;
      nacc++;
    }
    else
    {
      if(print_details)
      {
        Rcout << " - rejected." << std::endl;
      }
    }
  }
}



void MoveColonisationTime_NS(List &data, NumericVector parameters, double &loglik_cur, double w, IntegerVector Vq, IntegerVector Va, int &nacc_move)
{
  bool print_details = false;
  
  IntegerVector t_a = data["t_a"];
  IntegerVector t_d = data["t_d"];
  IntegerVector t_c = data["t_c"];
  IntegerVector source = data["source"];
  
  List data_can = clone(data);
  IntegerVector t_c_can = data_can["t_c"];
  IntegerVector source_can = data_can["source"];
  
  IntegerVector colonised_population = union_(Vq,Va);
  int target = SampleVector(colonised_population);
  //int last_day = ReturnLastDay(data, target);
  int last_day = t_d[target];
  
  IntegerMatrix screening_matrix = data["screening_matrix"];
  IntegerVector target_row = screening_matrix(target,_);
  IntegerVector positive_days = WhichVec(1, target_row);
  if(positive_days.length() > 0)
  {
    if(positive_days[0] < last_day)
    {
      last_day = positive_days[0];
    }
    
  }
  
  
  double log_prop_ratio = 0.0;
  double U = R::runif(0.0,1.0);
  if(U < w)
  {
    // Propose the target is an importation
    t_c_can[target] = t_a[target];
    source_can[target] = -1;
    if(source[target]==-1)
    {
      // Importation -> Importation
      // Do nothing, proposal ratio = 1
    }
    else
    {
      // Acquisition -> Importation
      log_prop_ratio = log(1-w) - log(w) - log(last_day-t_a[target]+1);
    }
  }
  else
  {
    // Propose the target is an acquisition
    //if(last_day < t_a[target]) return; // The target cannot be an acquisitioon
    t_c_can[target] = SampleVector(seq(t_a[target],last_day));
    IntegerVector possible_infectors = ReturnPossibleInfectors(data_can, target);
    if(possible_infectors.length()==0) return;
    source_can[target] = SampleVector(possible_infectors);
    if(source[target] != -1)
    {
      // Acquisition -> Acquisition
      //log_prop_ratio = log(total_col_pop(data_can, t_c_can[target])) - log(total_col_pop(data, t_c[target]));
    }
    else
    {
      // Importation -> Acquisition
      log_prop_ratio = log(w) + log(last_day - t_a[target] + 1) - log(1-w);
    }
  }
  
  //if(t_c[target] == t_c_can[target]) return;
  
  
  double loglik_can = CalculateLogLikelihood_NS(data_can, parameters);
  U = R::runif(0.0,1.0);
  //Rcout << "t_a = " << t_a[target] << ", t_d = " << t_d[target] << ", last_day = " << last_day << ", proposed coltime = " << t_c_can[target] << std::endl;
  
  if(print_details)
  {
    Rcout << std::endl;
    Rcout << "MOVING ID " << target << " from time " << t_c[target] << " to " << t_c_can[target] << " and from source " << source[target] << " to " << source_can[target] << std::endl;
    Rcout <<  "Loglik can = " << loglik_can << ", loglik cur = " << loglik_cur << ", log prop ratio = " << log_prop_ratio;
  }
  
  if(log(U) < loglik_can - loglik_cur + log_prop_ratio)
  {
    data = data_can;
    loglik_cur = loglik_can;
    nacc_move++;
    if(print_details) Rcout << " - accepted." << std::endl;
  }
  else
  {
    if(print_details) Rcout << " - rejected." << std::endl;
  }
}



void AddColonisationTime_NS(List &data, NumericVector parameters, double &loglik_cur, double w, IntegerVector Vs, IntegerVector &Va, int &nacc_add)
{
  bool print_details = false;
  IntegerVector t_a = data["t_a"];
  IntegerVector t_c = data["t_c"];
  IntegerVector t_d = data["t_d"];
  IntegerVector source = data["source"];
  
  List data_can = clone(data);
  IntegerVector t_c_can = data_can["t_c"];
  IntegerVector source_can = data_can["source"];
  
  IntegerVector susceptible_population = setdiff(Vs,Va);
  if(susceptible_population.length()==0) return;
  IntegerVector nadd0 = ReturnNadd0(data, Va);
  int target = SampleVector(susceptible_population);
  double log_prop_ratio = 0;
  
  double U = R::runif(0.0,1.0);
  if(U < w)
  {
    // Propose an importation
    t_c_can[target] = t_a[target];
    source_can[target] = -1;
    
    //log_prop_ratio = log(susceptible_population.length()) - log(w) - log(1+nadd0.length());
    log_prop_ratio = log(susceptible_population.length()) - log(w) - log(1+Va.length());
  }
  else
  {
    // Propose an acquisition
    t_c_can[target] = SampleVector(seq(t_a[target],t_d[target]));
    IntegerVector possible_infectors = ReturnPossibleInfectors(data_can, target);
    if(possible_infectors.length()==0) return;
    source_can[target] = SampleVector(possible_infectors);
    
    //log_prop_ratio = log(total_col_pop(data_can, t_c_can[target])) + log(susceptible_population.length()) + 
    //  log(t_d[target] - t_a[target] + 1) - log(1-w) - log(1+nadd0.length());
    log_prop_ratio = log(susceptible_population.length()) + log(t_d[target]-t_a[target]+1) - log(1-w) - log(1+Va.length());
  }
  
  //if(!CheckData(data_can)) stop("error with data with the adding step!");
  double loglik_can = CalculateLogLikelihood_NS(data_can, parameters);
  U = R::runif(0.0,1.0);
  
  
  if(print_details)
  {
    Rcout << std::endl;
    Rcout << "ADDING ID " << target << " at time " << t_c_can[target] << " with source " << source_can[target] << std::endl;
    Rcout << "Loglik can = " << loglik_can << ", loglik cur = " << loglik_cur << ", log prop ratio = " << log_prop_ratio;
  }
  if(log(U) < loglik_can - loglik_cur + log_prop_ratio)
  {
    data = data_can;
    loglik_cur = loglik_can;
    Va.push_back(target);
    nacc_add++;
    if(print_details) Rcout << " - accepted. Number added = " << Va.length() << std::endl;
  }
  else
  {
    if(print_details) Rcout << " - rejected." << std::endl;
  }
  
}

void RemoveColonisationTime_NS(List &data, NumericVector parameters, double &loglik_cur, double w, IntegerVector Vs, IntegerVector &Va, int &nacc_remove)
{
  bool print_details = false;
  IntegerVector t_a = data["t_a"];
  IntegerVector t_c = data["t_c"];
  IntegerVector t_d = data["t_d"];
  IntegerVector source = data["source"];
  
  List data_can = clone(data);
  IntegerVector t_c_can = data_can["t_c"];
  IntegerVector source_can = data_can["source"];
  
  IntegerVector susceptible_population = setdiff(Vs,Va);
  //IntegerVector nadd0 = ReturnNadd0(data, Va);
  //if(nadd0.length()==0) return;
  if(Va.length()==0) return;
  int target = SampleVector(Va);
  double log_prop_ratio = 0;
  
  source_can[target] = -2;
  t_c_can[target] = -1;
  
  if(source[target]==-1)
  {
    // Removing an importation
    //log_prop_ratio = log(nadd0.length()) + log(w) - log(susceptible_population.length()+1);
    log_prop_ratio = log(Va.length()) + log(w) - log(susceptible_population.length()+1);
  }
  else
  {
    // Removing an acquisition
    //log_prop_ratio = log(nadd0.length()) + log(1-w) - log(t_d[target]-t_a[target]+1) - log(susceptible_population.length()+1) - log(total_col_pop(data,t_c[target]));
    log_prop_ratio = log(Va.length()) + log(1-w) - log(t_d[target]-t_a[target]+1) - log(susceptible_population.length()+1);// - log(col_pop(data,t_c[target]));
    
  }
  
  //if(!CheckData(data_can)) stop("error with data in the remove step!");
  //Rcout << "Before likelihood" << std::endl;
  double loglik_can = CalculateLogLikelihood_NS(data_can, parameters);
  //Rcout << "After likelihood" << std::endl;
  double U = R::runif(0.0,1.0);
  if(print_details)
  {
    Rcout << std::endl;
    Rcout << "REMOVING ID " << target << " from time " << t_c[target] << " with source " << source[target] << std::endl;
    Rcout << "Loglik can = " << loglik_can << ", loglik cur = " << loglik_cur << ", log prop ratio = " << log_prop_ratio;
  }
  if(log(U) < loglik_can - loglik_cur + log_prop_ratio)
  {
    data = data_can;
    loglik_cur = loglik_can;
    int Va_loc = WhichVec(target, Va)[0];
    Va = RemoveElement(Va_loc, Va);
    if(print_details) Rcout << " - accepted. Number added = " << Va.length() << std::endl;
    nacc_remove++;
  }
  else
  {
    if(print_details) Rcout << " - rejected." << std::endl;
  }
}


// Update the imputed nodes by resimulating all of them  
List UpdateImputedNodesMoveAll(List data_cur, List data_can, NumericVector parameters)
{
  bool print_details = false;
  
  if(print_details) Rcout << "Attempt to return current and candidate nodes to impute" << std::endl;
  List nodes_to_impute_cur = ReturnNodesToImpute(data_cur);
  List nodes_to_impute_can = ReturnNodesToImpute(data_can);
  if(print_details) Rcout << "Returned nodes to impute" << std::endl;
  
  // check if the nodes are the same
  
  
  
  /*
   bool all_equal = false;
   if(imputed_genetic_ids_cur.length() == imputed_genetic_ids_can.length())
   {
   for(int i = 0; i<imputed_genetic_ids_can.length(); i++)
   {
   if((imputed_genetic_ids_cur[i] != imputed_genetic_ids_can[i]) || (imputed_sample_times_cur[i] != imputed_sample_times_can[i])
   || (imputed_variant_numbers_cur[i] != imputed_variant_numbers_can[i]))
   {
   // The nodes to impute are the same as before - now check they all have parents and children
   for(int j = 0; j<imputed_genetic_ids_cur.length(); j++)
   {
   int current_node = imputed_genetic_ids_cur[j];
   }
   
   
   all_equal = false;
   break;
   }
   }
   if(all_equal)
   {
   if(print_details)
   {
   Rcout << "Same nodes as before, do not update distances" << std::endl;
   }
   data_can["genetic_contribution"] = 0;
   return data_can;
   }
   }
   */
  
  
  
  
  double genetic_contribution = 0;
  
  // look at the reverse process and determine the genetic contribution from the imputed nodes
  IntegerVector imputed_nodes_cur = data_cur["imputed_nodes"];
  IntegerVector genetic_ids_cur = data_cur["genetic_ids"];
  IntegerVector sample_times_cur = data_cur["sample_times"];
  IntegerVector variant_numbers_cur = data_cur["variant_numbers"];
  IntegerVector t_c_cur = data_cur["t_c"];
  IntegerVector source_cur = data_cur["source"];
  IntegerMatrix genetic_matrix_cur = data_cur["genetic_matrix"];
  
  IntegerVector observed_idx = WhichVec(0,imputed_nodes_cur);
  IntegerVector observed_genetic_ids = as<IntegerVector>(genetic_ids_cur[observed_idx]);
  IntegerVector observed_sample_times = as<IntegerVector>(sample_times_cur[observed_idx]);
  IntegerVector observed_variant_numbers = as<IntegerVector>(variant_numbers_cur[observed_idx]);
  
  IntegerVector updated_genetic_ids = clone(observed_genetic_ids);
  IntegerVector updated_sample_times = clone(observed_sample_times);
  IntegerVector updated_variant_numbers = clone(observed_variant_numbers);
  
  
  // Check if the genetic source vectors are the same, if so do not update nodes
  IntegerVector imputed_genetic_ids_cur = nodes_to_impute_cur["nodes"];
  IntegerVector imputed_sample_times_cur = nodes_to_impute_cur["times"];
  IntegerVector imputed_variant_numbers_cur = nodes_to_impute_cur["variant_numbers"];
  
  IntegerVector imputed_genetic_ids_can = nodes_to_impute_can["nodes"];
  IntegerVector imputed_sample_times_can = nodes_to_impute_can["times"];
  IntegerVector imputed_variant_numbers_can = nodes_to_impute_can["variant_numbers"];
  
  
  // THIS NEEDS TO BE FIXED -- POTENTIALLY A PROBLEM WITH THE SAME GENETIC SOURCE VECTOR BUT DIFFERENT TIMES (?)
  /*
   if(imputed_genetic_ids_cur.length() == imputed_genetic_ids_can.length())
   {
   bool all_equal = true;
   IntegerVector updated_genetic_ids_check = clone(observed_genetic_ids);
   IntegerVector updated_sample_times_check = clone(observed_sample_times);
   IntegerVector updated_variant_numbers_check = clone(observed_variant_numbers);
   for(int i = 0; i<imputed_genetic_ids_can.length(); i++)
   {
   updated_genetic_ids_check.push_back(imputed_genetic_ids_can[i]);
   updated_sample_times_check.push_back(imputed_sample_times_can[i]);
   updated_variant_numbers_check.push_back(imputed_variant_numbers_can[i]);
   }
   IntegerVector t_c_can = data_can["t_c"];
   IntegerVector source_can = data_can["source"];
   
   List temp_data = List::create(Named("t_c") = t_c_can,
   Named("source") = source_can,
   Named("genetic_ids") = updated_genetic_ids_check,
   Named("sample_times") = updated_sample_times_check,
   Named("variant_numbers") = updated_variant_numbers_check);
   
   IntegerVector gen_source_cur = ReturnGenSourceVector_WHD(data_cur);
   IntegerVector gen_source_can = ReturnGenSourceVector_WHD(data_can);
   
   for(int i = 0; i<gen_source_cur.length(); i++)
   {
   if(gen_source_cur[i] != gen_source_can[i])
   {
   all_equal = false;
   break;
   }
   }
   if(all_equal)
   {
   data_can["genetic_contribution"] = 0;
   if(print_details)
   {
   double x = data_can["genetic_contribution"];
   Rcout << "Same nodes as before, do not update distances" << std::endl;
   Rcout << "Genetic contribution before = 0 is " <<  x << std::endl;
   }
   
   return data_can;
   }
   }
   */
  
  
  // print data to file for debug
  bool print_debug_file = false;
  if(print_debug_file)
  {
    // augmented data file
    std::string debug_file = "debug_file.dat";
    remove(debug_file.c_str());
    std::ofstream myfile6;
    myfile6.open(debug_file.c_str());
    assert(myfile6.is_open());
    
    
    PrintIntVectorToFile(myfile6, data_can["t_c"]);
    PrintIntVectorToFile(myfile6, data_can["source"]);
    PrintIntVectorToFile(myfile6, data_can["genetic_ids"]);
    PrintIntVectorToFile(myfile6, data_can["sample_times"]);
    PrintIntVectorToFile(myfile6, data_can["variant_numbers"]);
    PrintIntVectorToFile(myfile6, data_can["imputed_nodes"]);
    myfile6 << std::endl;
    myfile6.close();
  }
  
  
  
  if(print_details)
  {
    Rcout << "Attempt to calculate the contribution of the reverse process" << std::endl;
    Rcout << "imputed_genetic_ids_cur = " << imputed_genetic_ids_cur << std::endl;
    Rcout << "imputed_sample_times_cur = " << imputed_sample_times_cur << std::endl;
  }
  
  for(int i = 0; i<imputed_genetic_ids_cur.length(); i++)
  {
    int current_id = imputed_genetic_ids_cur[i];
    int current_time = imputed_sample_times_cur[i];
    int current_variant = imputed_variant_numbers_cur[i];
    
    updated_genetic_ids.push_back(current_id);
    updated_sample_times.push_back(current_time);
    updated_variant_numbers.push_back(current_variant);
    List temp_data = List::create(Named("t_c") = t_c_cur,
                                  Named("source") = source_cur,
                                  Named("genetic_ids") = updated_genetic_ids,
                                  Named("sample_times") = updated_sample_times,
                                  Named("variant_numbers") = updated_variant_numbers);
    
    if(print_details)
    {
      Rcout << "Current id = " << current_id << ", current time = " << current_time << ", current variant = " << current_variant << std::endl;
    }
    
    
    IntegerVector gen_source = ReturnGenSourceVector(temp_data);
    int imputed_node = gen_source.length()-1;
    //Rcout << "Gen source length = " << gen_source.length() << std::endl;
    //Rcout << "Gen source = " << gen_source << std::endl;
    //Rcout << "updated_genetic_ids = " << updated_genetic_ids << std::endl;
    //Rcout << "updated_sample_times = " << updated_sample_times << std::endl;
    int parent_node = gen_source[imputed_node];
    IntegerVector children_nodes = WhichVec(imputed_node, gen_source);
    
    if(parent_node == -1)
    {
      // the node is exterior
      int min_distance = 1000000000;
      int min_i = -1;
      int min_j = -1;
      for(int i = 0; i<(children_nodes.length()-1); i++)
      {
        for(int j = (i+1); j<children_nodes.length(); j++)
        {
          int node_i = children_nodes[i];
          int node_j = children_nodes[j];
          int current_distance = genetic_matrix_cur(node_i,node_j);
          if(current_distance < min_distance)
          {
            min_distance = current_distance;
            min_i = node_i;
            min_j = node_j;
          }
        }
      }
      
      if(min_i == -1 || min_j == -1)
      {
        Rcout << "Parent node = " << parent_node << ", imputed node = " << imputed_node << ", children nodes = " << children_nodes << std::endl;
        stop("Problem trying to find the minimum pairwise distance when removing nodes!");
      }
      
      //double upper = sample_times_cur[min_i] - sample_times_cur[imputed_node];
      //double lower = sample_times_cur[min_i] + sample_times_cur[min_j] - 2*sample_times_cur[imputed_node];
      double upper = updated_sample_times[min_i] - updated_sample_times[imputed_node];
      double lower = updated_sample_times[min_i] + updated_sample_times[min_j] - 2*updated_sample_times[imputed_node];
      double time_ratio = upper/lower;
      if(time_ratio < 0 || time_ratio > 1)
      {
        Rcout << "Time ratio = " << time_ratio << std::endl;
        stop("invalid time ratio");
      }
      
      int observed_distance = genetic_matrix_cur(imputed_node, min_i);
      if(print_details)
      {
        Rcout << "Removing exterior node = " << imputed_node << ", parent node = " << parent_node << ", children nodes = " 
              << children_nodes << ", min child = " << min_i << std::endl;
        Rcout << "Observed distance = " << observed_distance << ", min distance = " << min_distance << ", time ratio = " 
              << time_ratio << ", genetic contribution = "<< R::dbinom(observed_distance, min_distance, time_ratio, 1) << std::endl;
        Rcout << std::endl;
      }
      genetic_contribution += R::dbinom(observed_distance, min_distance, time_ratio, 1);
    }
    else
    {
      // interior node
      // We are removing an interior node
      int min_child = children_nodes[0];
      int min_distance = genetic_matrix_cur(parent_node, min_child);
      for(int j = 1; j<children_nodes.length(); j++)
      {
        int current_child = children_nodes[j];
        int current_distance = genetic_matrix_cur(parent_node, current_child);
        if(current_distance < min_distance)
        {
          min_child = current_child;
          min_distance = current_distance;
        }
      }
      double upper = updated_sample_times[imputed_node] - updated_sample_times[parent_node];
      double lower = updated_sample_times[min_child] - updated_sample_times[parent_node];
      double time_ratio = (double)(upper)/(double)(lower);
      int observed_distance = genetic_matrix_cur(parent_node, imputed_node);
      if(print_details)
      {
        Rcout << "Removing interior node = " << imputed_node << ", parent node = " << parent_node << ", children nodes = " 
              << children_nodes << ", min child = " << min_child << std::endl;
        Rcout << "Observed distance = " << observed_distance << ", min distance = " << min_distance << ", time ratio = " 
              << time_ratio << ", genetic contribution = " << R::dbinom(observed_distance, min_distance, time_ratio, 1) << std::endl;
        Rcout << std::endl;
      }
      /*
       if(print_details)
       {
       Rcout << "Observed distance = " << observed_distance << ", min distance = " << min_distance << ", time ratio = " << time_ratio << std::endl;
       Rcout << "sample_times_cur[imputed_node] = " << sample_times_cur[imputed_node] << ", sample_times_cur[parent_node] = " 
             << sample_times_cur[parent_node] << ", sample_times_cur[min_child] = " << sample_times_cur[min_child] << std::endl;
       }
       */
      genetic_contribution += R::dbinom(observed_distance, min_distance, time_ratio, 1);
    }
  }
  
  
  // Now simulate imputed nodes under the new configuration
  if(print_details)
  {
    Rcout << "Attempt to calculate the contribution of the forward process" << std::endl;
  }
  
  
  if(print_details)
  {
    Rcout << "Nodes cur = " << imputed_genetic_ids_cur << std::endl;
    Rcout << "Times cur = " << imputed_sample_times_cur << std::endl;
    Rcout << "Nodes can = " << imputed_genetic_ids_can << std::endl;
    Rcout << "Times can = " << imputed_sample_times_can << std::endl;
  }
  
  updated_genetic_ids = clone(observed_genetic_ids);
  updated_sample_times = clone(observed_sample_times);
  updated_variant_numbers = clone(observed_variant_numbers);
  
  
  
  int observed_length = observed_genetic_ids.length();
  IntegerMatrix observed_genetic_matrix(observed_length);
  for(int i = 0; i<observed_length; i++)
  {
    for(int j = 0; j<observed_length; j++)
    {
      observed_genetic_matrix(i,j) = genetic_matrix_cur(i,j);
    }
  }
  
  IntegerMatrix updated_genetic_matrix = clone(observed_genetic_matrix);
  IntegerVector t_c_can = data_can["t_c"];
  IntegerVector source_can = data_can["source"];
  
  //Rcout << std::endl;
  //PrintMatrix(updated_genetic_matrix);
  //Rcout << std::endl;
  
  for(int i = 0; i<imputed_genetic_ids_can.length(); i++)
  {
    int current_id = imputed_genetic_ids_can[i];
    int current_time = imputed_sample_times_can[i];
    int current_variant = imputed_variant_numbers_can[i];
    
    updated_genetic_ids.push_back(current_id);
    updated_sample_times.push_back(current_time);
    updated_variant_numbers.push_back(current_variant);
    List temp_data = List::create(Named("t_c") = t_c_can,
                                  Named("source") = source_can,
                                  Named("genetic_ids") = updated_genetic_ids,
                                  Named("sample_times") = updated_sample_times,
                                  Named("variant_numbers") = updated_variant_numbers);
    if(print_details)
    {
      Rcout << "Current id = " << current_id << ", current time = " << current_time << ", current variant = " << current_variant << std::endl;
    }
    IntegerVector gen_source = ReturnGenSourceVector(temp_data);
    IntegerMatrix genetic_matrix_can(gen_source.length());
    
    // copy matrix entries
    for(int i = 0; i<updated_genetic_matrix.nrow(); i++)
    {
      for(int j = 0; j<updated_genetic_matrix.ncol(); j++)
      {
        genetic_matrix_can(i,j) = updated_genetic_matrix(i,j);
      }
    }
    
    int imputed_node = gen_source.length()-1;
    int parent_node = gen_source[imputed_node];
    IntegerVector children_nodes = WhichVec(imputed_node, gen_source);
    if(parent_node == -1)
    {
      // impute an exterior node
      // find the minimum of the pairwise distances
      int min_distance = 1000000000;
      int min_i = -1;
      int min_j = -1;
      //Rcout << "Children nodes = " << children_nodes << std::endl;
      for(int i = 0; i<(children_nodes.length()-1); i++)
      {
        for(int j = (i+1); j<children_nodes.length(); j++)
        {
          int node_i = children_nodes[i];
          int node_j = children_nodes[j];
          int current_distance = updated_genetic_matrix(node_i,node_j);
          if(current_distance < min_distance)
          {
            min_distance = current_distance;
            min_i = node_i;
            min_j = node_j;
          }
        }
      }
      
      if(min_i == -1 || min_j == -1)
      {
        Rcout << std::endl << "**** Debug summary ****" << std::endl;
        Rcout << "Parent node = " << parent_node << ", imputed node = " << imputed_node << ", children nodes = " << children_nodes << std::endl;
        Rcout << "Gen source = " << gen_source << std::endl;
        Rcout << "Genetic ids = " << updated_genetic_ids << std::endl;
        Rcout << "Sample times = " << updated_sample_times << std::endl;
        IntegerVector colonised_individuals = WhichNotEqualVec(-1, t_c_can);
        Rcout << "Colonised individuals = " << colonised_individuals << std::endl;
        Rcout << "t_c_can = " << as<IntegerVector>(t_c_can[colonised_individuals]) << std::endl;
        Rcout << "source_can = " << as<IntegerVector>(source_can[colonised_individuals]) << std::endl;
        stop("Problem trying to find the minimum pairwise distance when adding nodes!");
      }
      
      
      // simulate the distance between min node i and the imputed node
      double upper = updated_sample_times[min_i] - updated_sample_times[imputed_node];
      double lower = updated_sample_times[min_i] + updated_sample_times[min_j] - 2*updated_sample_times[imputed_node];
      double time_ratio = upper/lower;
      if(time_ratio < 0 || time_ratio > 1)
      {
        Rcout << "Time ratio = " << time_ratio << std::endl;
        stop("invalid time ratio");
      }
      
      int simulated_distance = R::rbinom(min_distance, time_ratio);
      genetic_contribution -= R::dbinom(simulated_distance, min_distance, time_ratio, 1);
      if(print_details)
      {
        Rcout << "Adding exterior node = " << imputed_node << ", parent node = " << parent_node << ", children nodes = " 
              << children_nodes << ", min child = " << min_i << std::endl;
        Rcout << "Simulated distance = " << simulated_distance << ", min distance = " << min_distance << ", time ratio = " 
              << time_ratio << ", genetic contribution = "<< R::dbinom(simulated_distance, min_distance, time_ratio, 1) << std::endl;
        Rcout << std::endl;
      }
      //Rcout << "Min i = " << min_i << ", min j = " << min_j << ", simulated distance = " << simulated_distance << std::endl;  
      genetic_matrix_can(min_i, imputed_node) = simulated_distance;
      genetic_matrix_can(imputed_node, min_i) = simulated_distance;
      for(int i = 0; i<(gen_source.length()-1); i++)
      {
        IntegerVector loc_in_child_nodes = WhichVec(i, children_nodes);
        if(i != min_i)
        {
          // need to double check that these make sense
          if(loc_in_child_nodes.length()==0)
          {
            genetic_matrix_can(i, imputed_node) = genetic_matrix_can(i, min_i) + simulated_distance;
            genetic_matrix_can(imputed_node, i) = genetic_matrix_can(min_i, i) + simulated_distance;
          }
          else
          {
            genetic_matrix_can(i, imputed_node) = genetic_matrix_can(i, min_i) - simulated_distance;
            genetic_matrix_can(imputed_node, i) = genetic_matrix_can(min_i, i) - simulated_distance;
          }
        }
      }
    }
    else
    {
      // impute an interior node
      // Find the child that gives the minimum distance
      int min_child = children_nodes[0];
      int min_distance = genetic_matrix_can(parent_node, min_child);
      for(int i = 1; i<children_nodes.length(); i++)
      {
        int current_child = children_nodes[i];
        int current_distance = genetic_matrix_can(parent_node, current_child);
        if(current_distance < min_distance)
        {
          min_distance = current_distance;
          min_child = current_child;
        }
      }
      
      double upper = (double)(updated_sample_times[imputed_node]) - (double)(updated_sample_times[parent_node]);
      double lower = (double)(updated_sample_times[min_child]) - (double)(updated_sample_times[parent_node]);
      double time_ratio = upper/lower;
      
      if(time_ratio < 0 || time_ratio > 1)
      {
        Rcout << "Time ratio = " << time_ratio << std::endl;
        stop("invalid time ratio");
      }
      
      // Simulate a distance
      int simulated_distance = R::rbinom(min_distance, time_ratio);
      genetic_contribution -= R::dbinom(simulated_distance, min_distance, time_ratio, 1);
      if(print_details)
      {
        Rcout << "Adding interior node = " << imputed_node << ", parent node = " << parent_node << ", children nodes = " 
              << children_nodes << ", min child = " << min_child << std::endl;
        Rcout << "Simulated distance = " << simulated_distance << ", min distance = " << min_distance << ", time ratio = " 
              << time_ratio << ", genetic contribution = "<< R::dbinom(simulated_distance, min_distance, time_ratio, 1) << std::endl;
        Rcout << std::endl;
      }
      genetic_matrix_can(imputed_node, parent_node) = simulated_distance;
      genetic_matrix_can(parent_node, imputed_node) = simulated_distance;
      
      
      for(int i = 0; i<(gen_source.length()-1); i++)
      {
        if(i != parent_node)
        {
          IntegerVector children_loc = WhichVec(i,children_nodes);
          if(children_loc.length()==0)
          {
            // not a child node, calculate distance by summing 
            int dist = CalculateDistanceBetweenNodes_Rcpp(i, imputed_node, gen_source, genetic_matrix_can);
            genetic_matrix_can(i,imputed_node) = dist;
            genetic_matrix_can(imputed_node,i) = dist;
          }
          else
          {
            // child node, subtract
            genetic_matrix_can(i, imputed_node) = genetic_matrix_can(i, parent_node) - simulated_distance;
            genetic_matrix_can(imputed_node, i) = genetic_matrix_can(parent_node, i) - simulated_distance;
          }
        }
      }
      
      
      
      
      /*
       // fill in the rest 
       // ReturnPathToRoot(IntegerVector gen_source, int node)
       if(print_details) Rcout << "Gen source = " << gen_source << std::endl;
       if(print_details) Rcout << "Genetic ids =  " << updated_genetic_ids << std::endl;
       if(print_details) Rcout << "sample times = " << updated_sample_times << std::endl;
       for(int i = 0; i<(gen_source.length()-1); i++)
       {
       IntegerVector path_to_root = ReturnPathToRoot_Rcpp(gen_source, i);
       //Rcout << "Path to root = " << path_to_root << std::endl;
       
       IntegerVector imputed_loc_in_path = WhichVec(imputed_node, path_to_root);
       if(imputed_loc_in_path.length() == 0)
       {
       // The path does not contain the imputed node
       genetic_matrix_can(i, imputed_node) = genetic_matrix_can(i, parent_node) + simulated_distance;
       genetic_matrix_can(imputed_node, i) = genetic_matrix_can(parent_node, i) + simulated_distance;
       }
       else
       {
       // the path does contain the imputed node
       genetic_matrix_can(i, imputed_node) = genetic_matrix_can(i, parent_node) - simulated_distance;
       genetic_matrix_can(imputed_node, i) = genetic_matrix_can(parent_node, i) - simulated_distance;
       }
       }
       
       */
    }
    
    
    
    updated_genetic_matrix = genetic_matrix_can;
    //PrintMatrix(updated_genetic_matrix);
  }
  
  bool negative_matrix = DoesMatrixContainNegativeEntries(updated_genetic_matrix);
  if(negative_matrix)
  {
    //Rcout << "gen source can = " << gen_source_can << std::endl;
    //PrintMatrix(genetic_matrix_can);
    
    //Rcout << "Children nodes = " << children_nodes << std::endl;
    Rcout << "Negative moves in the matrix - no move is made!" << std::endl;
    genetic_contribution -= 1000; // make the move impossible
    //stop("negative values in the matrix!");
  }
  
  IntegerVector t_a = data_can["t_a"];
  IntegerVector t_d = data_can["t_d"];
  IntegerMatrix screening_matrix = data_can["screening_matrix"];
  int sequence_length = data_can["sequence_length"];
  
  
  IntegerVector updated_imputed_nodes(updated_genetic_ids.length());
  for(int i = observed_genetic_ids.length(); i<updated_genetic_ids.length(); i++)
  {
    updated_imputed_nodes[i] = 1;
  }
  if(print_details) Rcout << "Updated imputed nodes = " << updated_imputed_nodes << std::endl;
  
  
  // fix any imputed master distances that may be missing
  List temp_data = List::create(Named("t_c") = t_c_can,
                                Named("source") = source_can,
                                Named("genetic_ids") = updated_genetic_ids,
                                Named("sample_times") = updated_sample_times,
                                Named("variant_numbers") = updated_variant_numbers);
  
  IntegerVector gen_source_cur = ReturnGenSourceVector(data_cur);
  IntegerVector gen_source_can = ReturnGenSourceVector(temp_data);
  
  if(print_details) Rcout << "gen source cur = " << gen_source_cur << std::endl;
  if(print_details) Rcout << "gen source can = " << gen_source_can << std::endl;
  
  int num_patients = data_cur["num_patients"];
  IntegerVector hcw_ind = data_cur["hcw_ind"];
  
  
  //Rcout << "Need to impute" << std::endl;
  List final_data = List::create(Named("t_a") = t_a, 
                                 Named("t_c") = t_c_can,
                                 Named("t_d") = t_d,
                                 Named("source") = source_can,
                                 Named("screening_matrix") = screening_matrix,
                                 Named("genetic_ids") = updated_genetic_ids,
                                 Named("sample_times") = updated_sample_times,
                                 Named("variant_numbers") = updated_variant_numbers,
                                 Named("genetic_matrix") = updated_genetic_matrix,
                                 Named("sequence_length") = sequence_length,
                                 Named("imputed_nodes") = updated_imputed_nodes,
                                 Named("genetic_contribution") = genetic_contribution,
                                 Named("num_patients") = num_patients,
                                 Named("hcw_ind") = hcw_ind);
  
  return final_data;
}


/* Given an infection event, does a node need to be imputed? the way to check this is to see if there is eventually a sequence in the path of
 * the infector and also the infection
 * 
 */
// [[Rcpp::export]]
bool DoesNodeNeedToBeImputed(List data, int target, int current_variant)
{
  bool print_details = false;
  //if(target==205) print_details = true;
  IntegerVector t_c = data["t_c"];
  IntegerVector source = data["source"];
  IntegerVector genetic_ids = data["genetic_ids"];
  IntegerVector sample_times = data["sample_times"];
  IntegerVector variant_numbers = data["variant_numbers"];
  
  IntegerVector variant_idx = WhichVec(current_variant,variant_numbers);
  
  
  int col_time = t_c[target];
  int infector = source[target];
  
  if(print_details) Rcout << "Target = " << target << ", col time = " << col_time << ", infector = " << infector << std::endl;
  // First check the future events of the target
  IntegerVector target_genetic_idx = WhichVec(target, genetic_ids);
  IntegerVector target_variant_genetic_idx = intersect(target_genetic_idx, variant_idx);
  IntegerVector target_sample_times = as<IntegerVector>(sample_times[target_variant_genetic_idx]);
  
  // check time of colonisation
  IntegerVector col_time_idx = WhichVec(col_time, target_sample_times);
  if(col_time_idx.length() > 0) return false;
  
  
  bool forward_sequence_found_target = DoesTargetHaveSequenceGreaterThanTime(data, col_time, target, current_variant);
  
  if(!forward_sequence_found_target)
  {
    // Could not find a further sequence in the target, check any children and their children
    IntegerVector primary_targets = WhichVec(target, source);
    bool still_searching_secondary = true;
    while(still_searching_secondary)
    {
      if(primary_targets.length()==0)
      {
        // no one found in this branch, node does not need to be imputed
        return false;
      }
      IntegerVector secondary_targets;
      for(int i = 0; i<primary_targets.length(); i++)
      {
        int current_target = primary_targets[i];
        forward_sequence_found_target = DoesTargetHaveSequenceGreaterThanTime(data, col_time, current_target, current_variant);
        if(!forward_sequence_found_target)
        {
          // nothing found for this target, store his children to check later
          IntegerVector current_target_children = WhichVec(current_target, source);
          for(int j = 0; j<current_target_children.length(); j++)
          {
            secondary_targets.push_back(current_target_children[j]);
          }
        }
        else
        {
          // found a sequence, we can exit early
          still_searching_secondary = false;
          break;
        }
      }
      primary_targets = secondary_targets;
    }
  }
  
  
  if(!forward_sequence_found_target) return false;
  
  
  // now do the same but with the infector
  bool forward_sequence_found_infector = false;
  
  // check if the infector has a swab at the time of infection, OR another infection with a swab at that time
  IntegerVector infector_genetic_idx = WhichVec(infector, genetic_ids);
  IntegerVector infector_variant_genetic_idx = intersect(infector_genetic_idx, variant_idx);
  IntegerVector infector_sample_times = as<IntegerVector>(sample_times[infector_variant_genetic_idx]);
  
  // check time of colonisation
  col_time_idx = WhichVec(col_time, infector_sample_times);
  if(col_time_idx.length() > 0) return false;
  
  // check sources time of col
  IntegerVector infector_sources = WhichVec(infector, source);
  
  if(print_details) Rcout << "infector sources = " << infector_sources << std::endl;
  IntegerVector primary_targets;
  for(int i = 0; i<infector_sources.length(); i++)
  {
    int current_source = infector_sources[i];
    if(print_details) Rcout << "i = " << i << ", current source = " << current_source << std::endl;
    if(current_source != target)
    {
      int current_source_col_time = t_c[current_source];
      // only look at individuals who are colonised at the same time or after the target
      if(current_source_col_time >= col_time)
      {
        if(current_source_col_time == col_time)
        {
          IntegerVector infector_other_genetic_idx = WhichVec(current_source, genetic_ids);
          IntegerVector infector_variant_other_genetic_idx = intersect(infector_other_genetic_idx, variant_idx);
          IntegerVector infector_other_sample_times = as<IntegerVector>(sample_times[infector_variant_other_genetic_idx]);
          if(print_details) Rcout << "infector_other_genetic_idx = " << infector_other_genetic_idx << ", infector_other_sample_times = " << infector_other_sample_times 
                                  << std::endl;
          col_time_idx = WhichVec(col_time, infector_other_sample_times);
          if(print_details) Rcout << "col time = " << col_time << std::endl;
          if(col_time_idx.length() > 0)
          {
            return false;
          }
        }
        primary_targets.push_back(current_source);
      }
      
    }
  }
  
  forward_sequence_found_infector = DoesTargetHaveSequenceGreaterThanTime(data, col_time, infector, current_variant);
  
  if(print_details) Rcout << "Infector = " << infector << ", forward seq found = " << forward_sequence_found_infector << std::endl;
  
  // look through children
  bool still_searching_secondary = true;
  if(forward_sequence_found_infector)
  {
    still_searching_secondary = false;
  }
  while(still_searching_secondary)
  {
    if(primary_targets.length()==0)
    {
      // no one found in this branch, node does not need to be imputed
      return false;
    }
    IntegerVector secondary_targets;
    if(print_details) Rcout << "Primary targets = " << primary_targets << std::endl;
    for(int i = 0; i<primary_targets.length(); i++)
    {
      int current_target = primary_targets[i];
      forward_sequence_found_infector = DoesTargetHaveSequenceGreaterThanTime(data, col_time, current_target, current_variant);
      if(print_details) Rcout << "Current target = " << current_target << ", foward sequence found = " << forward_sequence_found_infector << std::endl;
      if(!forward_sequence_found_infector)
      {
        // nothing found for this target, store his children to check later
        IntegerVector current_target_children = WhichVec(current_target, source);
        for(int j = 0; j<current_target_children.length(); j++)
        {
          secondary_targets.push_back(current_target_children[j]);
        }
      }
      else
      {
        // found a sequence, we can exit early
        still_searching_secondary = false;
        break;
      }
    }
    primary_targets = secondary_targets;
  }
  
  if(print_details)
  {
    Rcout << "Forward sequence found in infector = " << forward_sequence_found_infector << std::endl;
    Rcout << "Forward sequence found in target = " << forward_sequence_found_target << std::endl;
  }
  
  
  if(forward_sequence_found_infector && forward_sequence_found_target)
  {
    return true;
  }
  return false;
}

// [[Rcpp::export]]
List Rcpp_sort_imputed_nodes(IntegerVector nodes, IntegerVector times, IntegerVector variants)
{
  IntegerVector idx = seq_along(times) - 1;
  
  std::sort(idx.begin(), idx.end(), [&](int i, int j){
    if(nodes[i]==nodes[j]) {
      return variants[i] < variants[j];
    }
    return nodes[i] < nodes[j];
  });
  
  List return_data = List::create(Named("nodes") = nodes[idx],
                                  Named("times") = times[idx],
                                                        Named("variant_numbers") = variants[idx]);
  return return_data;
}


// [[Rcpp::export]]
List ReturnNodesToImpute(List data)
{
  bool print_details = false;
  IntegerVector t_c = data["t_c"];
  IntegerVector source = data["source"];
  IntegerVector genetic_ids = data["genetic_ids"];
  IntegerVector sample_times = data["sample_times"];
  IntegerVector variant_numbers = data["variant_numbers"];
  IntegerVector colonised_idx = WhichNotEqualVec(-1, t_c);
  
  //Timer timer;
  //timer.step("start");        // record the starting point
  
  
  //Rcout << "Returning nodes to impute " << std::endl;
  IntegerVector imputed_nodes(genetic_ids.length());
  if(data.containsElementNamed("imputed_nodes"))
  {
    imputed_nodes = data["imputed_nodes"];
  }
  
  // store the ID of all patients who were infected/acquisitions
  IntegerVector infection_events;
  for(int i = 0; i<colonised_idx.length(); i++)
  {
    if(source[colonised_idx[i]] != -1)
    {
      infection_events.push_back(colonised_idx[i]);
    }
  }
  
  //timer.step("Store infection events");      // record the first step
  
  IntegerVector infection_event_times = as<IntegerVector>(t_c[infection_events]);
  
  //Rcout << "Before sort - infection events = " << infection_events << ", times = " << infection_event_times << std::endl;
  
  // sort the events by times
  infection_events = Rcpp_sort(infection_events, infection_event_times);
  
  //timer.step("Sort infection events");      // record the next step
  //Rcout << "Infection events = " << infection_events << std::endl;
  //Rcout << "After sort - infection events = " << infection_events << ", times = " << infection_event_times << std::endl;
  IntegerVector updated_genetic_ids = clone(genetic_ids);
  IntegerVector updated_sample_times = clone(sample_times);
  IntegerVector updated_variant_numbers = clone(variant_numbers);
  IntegerVector imputed_idx = WhichVec(1, imputed_nodes);
  for(int i = 0; i<imputed_idx.length(); i++)
  {
    int current_idx = imputed_idx[i];
    updated_genetic_ids = RemoveElement(current_idx-i,updated_genetic_ids);
    updated_sample_times = RemoveElement(current_idx-i,updated_sample_times);
    updated_variant_numbers = RemoveElement(current_idx-i,updated_variant_numbers);
  }
  List temp_data = List::create(Named("t_c") = t_c,
                                Named("source") = source,
                                Named("genetic_ids") = updated_genetic_ids,
                                Named("sample_times") = updated_sample_times,
                                Named("variant_numbers") = updated_variant_numbers);
  IntegerVector nodes_to_impute;
  IntegerVector imputed_times;
  IntegerVector variant_numbers_to_impute;
  // go through each event and add a node if needed
  IntegerVector unique_variants = unique(variant_numbers);
  IntegerVector variants_to_search;
  for(int i = 0; i<unique_variants.length(); i++)
  {
    int current_variant = unique_variants[i];
    IntegerVector variant_idx = WhichVec(current_variant, variant_numbers);
    if(variant_idx.length() > 1)
    {
      variants_to_search.push_back(current_variant);
    }
  }
  //Rcout << "Variants to search length = " << variants_to_search.length() << ", unique variants length = " << unique_variants.length() << std::endl;
  
  for(int i = 0; i<infection_events.length(); i++)
  {
    for(int j = 0; j<variants_to_search.length(); j++)
    {
      int current_variant = variants_to_search[j];
      int current_infection_event = infection_events[i];
      //Rcout << "Current infection event = " << current_infection_event << std::endl;
      bool does_node_need_impute = DoesNodeNeedToBeImputed(temp_data, current_infection_event, current_variant);
      if(print_details)
      {
        Rcout << "Trying ID = " << current_infection_event << " at time " << t_c[current_infection_event] << " from source = " 
              << source[current_infection_event] << " with variant =  " << current_variant << ", node needed = " << does_node_need_impute << std::endl;
      }
      if(does_node_need_impute)
      {
        // Check for duplicate
        IntegerVector id_locs = WhichVec(source[current_infection_event], nodes_to_impute);
        IntegerVector time_locs = WhichVec(t_c[current_infection_event], imputed_times);
        IntegerVector variant_locs = WhichVec(current_variant, variant_numbers_to_impute);
        IntegerVector id_time_intersect = intersect(id_locs, time_locs);
        IntegerVector id_time_variant_intersect = intersect(id_time_intersect, variant_locs);
        if(print_details)
        {
          Rcout << "id_locs = " << id_locs << std::endl;
          Rcout << "time_locs = " << time_locs << std::endl;
          Rcout << "variant_locs = " << variant_locs << std::endl;
          Rcout << "id_time_variant_intersect = " << id_time_variant_intersect << std::endl;
        }
        
        if(id_time_variant_intersect.length()==0)
        {
          updated_genetic_ids.push_back(source[current_infection_event]);
          updated_sample_times.push_back(t_c[current_infection_event]);
          updated_variant_numbers.push_back(current_variant);
          nodes_to_impute.push_back(source[current_infection_event]);
          imputed_times.push_back(t_c[current_infection_event]);
          variant_numbers_to_impute.push_back(current_variant);
          temp_data = List::create(Named("t_c") = t_c,
                                   Named("source") = source,
                                   Named("genetic_ids") = updated_genetic_ids,
                                   Named("sample_times") = updated_sample_times,
                                   Named("variant_numbers") = updated_variant_numbers);
        }
        
        
      }
    }
    
    
  }
  //timer.step("Determine if nodes need imputing");      // record the next step
  
  //NumericVector res(timer);   // 
  //Rcout << res << std::endl;
  
  
  if(print_details)
  {
    Rcout << "nodes_to_impute = " << nodes_to_impute << std::endl;
    Rcout << "imputed_times = " << imputed_times << std::endl;
    Rcout << "variant_numbers_to_impute = " << variant_numbers_to_impute << std::endl;
  }
  
  
  List data_out = Rcpp_sort_imputed_nodes(nodes_to_impute, imputed_times, variant_numbers_to_impute);
  
  //List data_out = UniqueSort(nodes_to_impute, imputed_times, variant_numbers_to_impute);
  IntegerVector final_nodes = data_out["nodes"];
  IntegerVector final_times = data_out["times"];
  IntegerVector final_variants = data_out["variant_numbers"];
  
  
  List return_data = List::create(Named("nodes") = final_nodes,
                                  Named("times") = final_times,
                                  Named("variant_numbers") = final_variants);
  /*
   List return_data = List::create(Named("nodes") = nodes_to_impute,
   Named("times") = imputed_times,
   Named("variant_numbers") = variant_numbers_to_impute);
   */
  return return_data;
}



// [[Rcpp::export]]
void MCMC_EPI_NS(List MCMC_options, IntegerVector t_a, IntegerVector t_c, IntegerVector t_d, 
                 IntegerVector source, IntegerVector hcw_ind, IntegerMatrix screening_matrix)
{
  // Initialise data
  List data = InitialiseData(t_a, t_c, t_d, source, hcw_ind, screening_matrix);
  int num_patients = data["num_patients"];
  
  // Load MCMC options
  const int iterations = MCMC_options["iterations"];
  const int num_updates = MCMC_options["num_updates"];
  NumericVector prior_parameters = MCMC_options["prior_parameters"];
  NumericVector initial_chain_state = MCMC_options["initial_chain_state"];
  NumericVector parameters = clone(initial_chain_state);
  NumericVector proposal_variance = MCMC_options["proposal_variance"];
  IntegerVector debug_flags = MCMC_options["debug_flags"];
  
  std::string output_file = MCMC_options["output_file"];
  std::string source_file = MCMC_options["source_file"];
  std::string coltime_file = MCMC_options["coltime_file"];
  
  // Check if there are healthcare workers
  bool include_hcw = true;
  if(t_c.length() == num_patients)
  {
    include_hcw = false;
    parameters[3] = 0.0;
    Rcout << "No healthcare workers present - use standard model" << std::endl;
  }
  else
  {
    Rcout << "Healthcare workers present, include beta_h term" << std::endl;
  }
  
  
  // Calculate other variables
  int true_positives = CalculateTruePositives(data);
  int nacc_beta_p = 0;
  int nacc_beta_h = 0;
  int nacc_add = 0;
  int nacc_move = 0;
  int nacc_remove = 0;
  double loglik = 0.0;
  IntegerVector Va;
  IntegerVector Vs = CalculateVs(data);
  IntegerVector Vq = CalculateVq(data);
  IntegerVector augmented_moves_proposed(5);
  
  // Write to files
  // Output file
  remove(output_file.c_str());
  std::ofstream myfile; // Define output stream
  myfile.open(output_file.c_str()); // Open file
  assert(myfile.is_open());
  PrintNumVectorToFile(myfile, parameters);
  
  // Col time file
  remove(coltime_file.c_str());
  std::ofstream myfile3;
  myfile3.open(coltime_file.c_str());
  assert(myfile3.is_open());
  PrintColonisationTimeSumToFile(myfile3, data);
  //PrintIntVectorToFile(myfile3, data["t_c"]);
  
  
  
  
  // Begin MCMC
  Rcout << "Begin MCMC" << std::endl;
  for(int i = 1; i < iterations; i++)
  {
    // Update z by gibbs
    if(debug_flags[0]==0)
    {
      parameters[0] = R::rbeta(prior_parameters[0] + true_positives, prior_parameters[1] + CalculateFalseNegatives(data));
    }
    
    // Update p by gibbs
    if(debug_flags[1]==0)
    {
      int import_sum = WhichVec(-1, data["source"]).length();
      parameters[1] = R::rbeta(prior_parameters[2] + import_sum, prior_parameters[3] + num_patients - import_sum);
    }
    
    
    loglik = CalculateLogLikelihood_NS(data, parameters);
    
    // Update beta_P by metropolis hastings
    if(debug_flags[2]==0)
    {
      UpdateTransmissionRate_P_NS(data, parameters, loglik, proposal_variance[0], prior_parameters[4], nacc_beta_p);
    }
    
    // update beta_H by metropolis hastings
    if(debug_flags[3]==0 && include_hcw)
    {
      UpdateTransmissionRate_H_NS(data, parameters, loglik, proposal_variance[1], prior_parameters[5], nacc_beta_h);
    }
    
    // Augmented data updates
    
    // Augmented data updates
    if(debug_flags[4]==0)
    {
      for(int update_counter = 0; update_counter < num_updates; update_counter++)
      {
        double w = 0.3;
        //double w = 0.00001;
        double move = R::runif(0.0,3.0);
        augmented_moves_proposed[floor(move)]++;
        if(floor(move) < 1)
        {
          // Move a colonisation time
          //Rcout << "Attempt to move a time" << std::endl;
          MoveColonisationTime_NS(data, parameters, loglik, w, Vq, Va, nacc_move);
        }
        else if(floor(move) < 2)
        {
          // Add a colonisation time
          AddColonisationTime_NS(data, parameters, loglik, w, Vs, Va, nacc_add);
        }
        else if(floor(move) < 3)
        {
          // Remove a colonisation time
          RemoveColonisationTime_NS(data, parameters, loglik, w, Vs, Va, nacc_remove);
        }
        
        
        
        
      }
      // Write colonisation time sum to file
      //PrintIntVectorToFile(myfile2, data["source"]);
      PrintColonisationTimeSumToFile(myfile3, data);
      //PrintIntVectorToFile(myfile3, data["t_c"]);
      
    }
    
    
    
    
    // Write parameters to file
    PrintNumVectorToFile(myfile, parameters);
    
    
    if(i%100==0)
    {
      Rcout << "Iteraion " << i << " completed, number added by algorithm = " << Va.length() << std::endl;
    }
    
  }
  
  // Close files
  myfile.close();
  //myfile2.close();
  myfile3.close();
  
  double beta_h_prob = (double)nacc_beta_h/(double)iterations;
  double beta_p_prob = (double)nacc_beta_p/(double)iterations;
  double move_prob = (double)nacc_move/(double)(augmented_moves_proposed[0]);
  double add_prob = (double)nacc_add/(double)(augmented_moves_proposed[1]);
  double remove_prob = (double)nacc_remove/(double)(augmented_moves_proposed[2]);
  Rcout << "beta_P acceptance = " << beta_p_prob << ", beta_H acceptance = " << beta_h_prob << std::endl;
  Rcout << "Move prob = " << move_prob << ", add prob = " << add_prob << ", remove prob = " << remove_prob << std::endl;
  //Rcout << "Move prob = " << move_prob << ", add prob = " << add_prob << ", remove prob = " << remove_prob << std::endl;
}

List InitialiseMCMCImputedNodes2(List data)
{
  List nodes_to_impute = ReturnNodesToImpute(data);
  Rcout << "Returned nodes to impute" << std::endl;
  
  // look at the reverse process and determine the genetic contribution from the imputed nodes
  IntegerVector genetic_ids = data["genetic_ids"];
  IntegerVector sample_times = data["sample_times"];
  IntegerVector variant_numbers = data["variant_numbers"];
  IntegerVector t_c = data["t_c"];
  IntegerVector source = data["source"];
  IntegerMatrix genetic_matrix = data["genetic_matrix"];
  
  int original_length = genetic_ids.length();
  
  IntegerVector updated_genetic_ids = clone(genetic_ids);
  IntegerVector updated_sample_times = clone(sample_times);
  IntegerVector updated_variant_numbers = clone(variant_numbers);
  
  IntegerVector imputed_genetic_ids = nodes_to_impute["nodes"];
  IntegerVector imputed_sample_times = nodes_to_impute["times"];
  IntegerVector imputed_variant_numbers = nodes_to_impute["variant_numbers"];
  
  Rcout << "Nodes to impute = " << imputed_genetic_ids << std::endl;
  Rcout << "sample times to impute = " << imputed_sample_times << std::endl;
  Rcout << "variant numbers to impute = " << imputed_variant_numbers << std::endl;
  
  for(int i = 0; i<imputed_genetic_ids.length(); i++)
  {
    int current_id = imputed_genetic_ids[i];
    int current_time = imputed_sample_times[i];
    int current_variant = imputed_variant_numbers[i];
    
    genetic_ids.push_back(current_id);
    sample_times.push_back(current_time);
    variant_numbers.push_back(current_variant);
    
    List temp_data = List::create(Named("t_c") = t_c,
                                  Named("source") = source,
                                  Named("genetic_ids") = genetic_ids,
                                  Named("sample_times") = sample_times,
                                  Named("variant_numbers") = variant_numbers);
    IntegerVector gen_source = ReturnGenSourceVector(temp_data);
    int imputed_node = gen_source.length() - 1;
    int parent_node = gen_source[imputed_node];
    int node_to_copy = parent_node;
    if(parent_node == -1)
    {
      // imported case, copy first child
      IntegerVector children_nodes = WhichVec(imputed_node, gen_source);
      
      int min_distance = 1000000000;
      int min_i = -1;
      int min_j = -1;
      for(int i = 0; i<(children_nodes.length()-1); i++)
      {
        for(int j = (i+1); j<children_nodes.length(); j++)
        {
          int node_i = children_nodes[i];
          int node_j = children_nodes[j];
          int current_distance = genetic_matrix(node_i,node_j);
          if(current_distance < min_distance)
          {
            min_distance = current_distance;
            min_i = node_i;
            min_j = node_j;
          }
        }
      }
      
      if(min_i == -1 || min_j == -1)
      {
        Rcout << "Parent node = " << parent_node << ", imputed node = " << imputed_node << ", children nodes = " << children_nodes << std::endl;
        stop("Problem trying to find the minimum pairwise distance when removing nodes!");
      }
      
      node_to_copy = min_i;
    }
    
    IntegerMatrix updated_genetic_matrix(gen_source.length());
    for(int ii = 0; ii<(gen_source.length()-1); ii++)
    {
      for(int jj = 0; jj<(gen_source.length()-1); jj++)
      {
        updated_genetic_matrix(ii,jj) = genetic_matrix(ii,jj);
      }
      updated_genetic_matrix(ii,imputed_node) = genetic_matrix(ii,node_to_copy);
      updated_genetic_matrix(imputed_node,ii) = genetic_matrix(node_to_copy,ii);
    }
    
    genetic_matrix = updated_genetic_matrix;
  }
  
  IntegerVector imputed_nodes(genetic_ids.length());
  for(int i = original_length; i<genetic_ids.length(); i++)
  {
    imputed_nodes[i] = 1;
  }
  
  
  IntegerVector t_a = data["t_a"];
  IntegerVector t_d = data["t_d"];
  IntegerMatrix screening_matrix = data["screening_matrix"];
  int sequence_length = data["sequence_length"];
  int num_patients = data["num_patients"];
  IntegerVector hcw_ind = data["hcw_ind"];
  
  
  List final_data = List::create(Named("t_a") = t_a , 
                                 Named("t_c") = t_c,
                                 Named("t_d") = t_d,
                                 Named("source") = source,
                                 Named("screening_matrix") = screening_matrix,
                                 Named("genetic_ids") = genetic_ids,
                                 Named("sample_times") = sample_times,
                                 Named("variant_numbers") = variant_numbers,
                                 Named("genetic_matrix") = genetic_matrix,
                                 Named("sequence_length") = sequence_length,
                                 Named("imputed_nodes") = imputed_nodes,
                                 Named("num_patients") = num_patients,
                                 Named("hcw_ind") = hcw_ind);
  Rcout << "Before gen source" << std::endl;
  IntegerVector new_gen_source = ReturnGenSourceVector(final_data);
  Rcout << "After gen source " << std::endl;
  
  
  
  
  
  bool print_debug = true;
  if(print_debug)
  {
    Rcout << "Final gen source = " << new_gen_source << std::endl;
  }
  //stop("stop early");
  
  return final_data;
}


// [[Rcpp::export]]
void MCMC_EPI_GEN_HCW(List MCMC_options, int sequence_length, IntegerVector t_a, IntegerVector t_c, IntegerVector t_d, IntegerVector source, 
                        IntegerVector hcw_ind, IntegerMatrix screening_matrix, IntegerVector genetic_ids, IntegerVector sample_times,
                        IntegerVector variant_numbers, IntegerMatrix genetic_matrix)
{
  // Initialise data

  List data = InitialiseData_GEN(sequence_length, t_a, t_c, t_d, source, hcw_ind, screening_matrix, genetic_ids, sample_times, 
                                 variant_numbers, genetic_matrix);
  Rcout << "Data initialised" << std::endl;
  // Initialise Imputed nodes
  data = InitialiseMCMCImputedNodes2(data);
  

  
  
  int num_patients = data["num_patients"];
  
  // Load MCMC options
  const int iterations = MCMC_options["iterations"];
  const int num_updates = MCMC_options["num_updates"];
  NumericVector prior_parameters = MCMC_options["prior_parameters"];
  NumericVector initial_chain_state = MCMC_options["initial_chain_state"];
  NumericVector parameters = clone(initial_chain_state);
  NumericVector proposal_variance = MCMC_options["proposal_variance"];
  IntegerVector debug_flags = MCMC_options["debug_flags"];
  std::string output_file = MCMC_options["output_file"];
  std::string source_file = MCMC_options["source_file"];
  std::string coltime_file = MCMC_options["coltime_file"];
  std::string augmented_moves_file = "aug_moves.dat";
  std::string gen_source_file = MCMC_options["gen_source_file"];
  std::string sample_times_file = MCMC_options["sample_times_file"];
  std::string variant_numbers_file = MCMC_options["variant_numbers_file"];
  
  Rcout << "MCMC Options loaded" << std::endl;
  Rcout << "Initial chain state = " << initial_chain_state << std::endl;
  
  // Check if there are healthcare workers
  bool include_hcw = true;
  if(t_c.length() == num_patients)
  {
    include_hcw = false;
    parameters[3] = 0.0;
    Rcout << "No healthcare workers present - use standard model" << std::endl;
  }
  else
  {
    Rcout << "Healthcare workers present, include beta_h term" << std::endl;
  }
  
  // Calculate other variables
  int true_positives = CalculateTruePositives(data);
  int nacc_beta_p = 0;
  int nacc_beta_h = 0;
  int nacc_lambda = 0;
  int nacc_add = 0;
  int nacc_move = 0;
  int nacc_remove = 0;
  int nacc_change = 0;
  int nacc_swap = 0;
  //int nacc_move_multi = 0;
  int nacc_change_all = 0;
  int nacc_change_random = 0;
  double loglik = 0.0;
  IntegerVector Va;
  IntegerVector Vs = CalculateVs(data);
  IntegerVector Vq = CalculateVq(data);
  IntegerVector augmented_moves_proposed(7);
  // Write to files
  // Output file
  remove(output_file.c_str());
  std::ofstream myfile; // Define output stream
  myfile.open(output_file.c_str()); // Open file
  assert(myfile.is_open());
  PrintNumVectorToFile(myfile, parameters);
  
  // Source file
  remove(source_file.c_str());
  std::ofstream myfile2;
  myfile2.open(source_file.c_str());
  assert(myfile2.is_open());
  PrintIntVectorToFile(myfile2, data["source"]);
  
  // Col time file
  remove(coltime_file.c_str());
  std::ofstream myfile3;
  myfile3.open(coltime_file.c_str());
  assert(myfile3.is_open());
  PrintColonisationTimeSumToFile(myfile3, data);
  //PrintIntVectorToFile(myfile3, data["t_c"]);
  
  // augmented moves file
  remove(augmented_moves_file.c_str());
  std::ofstream myfile4;
  myfile4.open(augmented_moves_file.c_str());
  assert(myfile4.is_open());
  //PrintColonisationTimeSumToFile(myfile3, data);
  
  remove(gen_source_file.c_str());
  std::ofstream myfile5;
  myfile5.open(gen_source_file.c_str());
  assert(myfile5.is_open());
  IntegerVector gen_source = ReturnGenSourceVector(data);
  PrintIntVectorToFile(myfile5, gen_source);
  
  // Sample times file
  remove(sample_times_file.c_str());
  std::ofstream myfile6;
  myfile6.open(sample_times_file.c_str());
  assert(myfile6.is_open());
  PrintIntVectorToFile(myfile6, data["sample_times"]);
  
  // Variant numbers file
  remove(variant_numbers_file.c_str());
  std::ofstream myfile7;
  myfile7.open(variant_numbers_file.c_str());
  assert(myfile7.is_open());
  PrintIntVectorToFile(myfile7, data["variant_numbers"]);
  
  Rcout << "Output files initialised" << std::endl;
  // Begin MCMC
  Rcout << "Begin MCMC" << std::endl;
  for(int i = 1; i < iterations; i++)
  {
    // Update z by gibbs
    if(debug_flags[0]==0)
    {
      parameters[0] = R::rbeta(prior_parameters[0] + true_positives, prior_parameters[1] + CalculateFalseNegatives(data));
    }
    
    // Update p by gibbs
    if(debug_flags[1]==0)
    {
      int import_sum = CalculateImportationSum(data);
      parameters[1] = R::rbeta(prior_parameters[2] + import_sum, prior_parameters[3] + num_patients - import_sum);
    }
    
    List import_distance_info = CalculateImportDistance(data, parameters);
    double distance_sum = import_distance_info["distance_sum"];
    double import_counter = import_distance_info["import_counter"];
    
    if(debug_flags[5]==0)
    {
      NumericVector gamma_draws = Rcpp::rgamma(1, prior_parameters[7] + distance_sum, 1/(prior_parameters[8] + import_counter));
      parameters[5] = gamma_draws[0];
    }
    
    
    loglik = CalculateLogLikelihood(data, parameters);
    //Rcout << "Log likelihood calculated" << std::endl;
    
    
    // Update beta_P by metropolis hastings
    if(debug_flags[2]==0)
    {
      UpdateTransmissionRate_P(data, parameters, loglik, proposal_variance[0], prior_parameters[4], nacc_beta_p);
    }
    
    //Rcout << "Beta P updated" << std::endl;
    
    // update beta_H by metropolis hastings
    if(debug_flags[3]==0 && include_hcw)
    {
      UpdateTransmissionRate_H(data, parameters, loglik, proposal_variance[1], prior_parameters[5], nacc_beta_h);
    }
    //Rcout << "Beta H updated" << std::endl;
    
    if(debug_flags[4]==0)
    {
      UpdateMutationRate(data, parameters, loglik, proposal_variance[2], prior_parameters[6], nacc_lambda);
    }
    //Rcout << "Lambda updated" << std::endl;
    
    

    //Rcout << "mu updated" << std::endl;

    
    // Augmented data updates
    if(debug_flags[6]==0)
    {
      for(int update_counter = 0; update_counter < num_updates; update_counter++)
      {
        double w = 0.3;
        //double w = 0.00001;
        double move = R::runif(0.0,6.0);
        augmented_moves_proposed[floor(move)]++;
        if(floor(move) < 1)
        {
          // Move a colonisation time
          data = MoveColonisationTimeGenetic(data, parameters, loglik, w, Vq, Va, nacc_move, myfile4);
        }
        else if(floor(move) < 2)
        {
          // Add a colonisation time
          AddColonisationTime(data, parameters, loglik, w, Vs, Va, nacc_add);
        }
        else if(floor(move) < 3)
        {
          // Remove a colonisation time
          RemoveColonisationTime(data, parameters, loglik, w, Vs, Va, nacc_remove);
        }
        else if(floor(move) < 4)
        {
          // Change a source without changing the time of colonisation
          data = UpdateSource(data, parameters, loglik, w, Vq, Va, nacc_change, myfile4);
        }
        else if(floor(move) < 5)
        {
          // Swap a target with their source
          data = SwapTargetOffspring(data, parameters, loglik, w, Vq, Va, nacc_swap, myfile4);          
        }
        else if(floor(move) < 6)
        {
          data = UpdateRandomSources(data, parameters, loglik, w, Vq, Va, nacc_change_random, myfile4);
          
        }
        else if(floor(move) < 7)
        {
          data = UpdateAllSources(data, parameters, loglik, w, Vq, Va, nacc_change_all, myfile4);
        }
        
        
        
        
      }
      // Write colonisation time sum to file
      PrintIntVectorToFile(myfile2, data["source"]);
      PrintColonisationTimeSumToFile(myfile3, data);
      gen_source = ReturnGenSourceVector(data);
      PrintIntVectorToFile(myfile5, gen_source);
      PrintIntVectorToFile(myfile6, data["sample_times"]);
      PrintIntVectorToFile(myfile7, data["variant_numbers"]);
      //PrintIntVectorToFile(myfile3, data["t_c"]);
      
    }
    
    
    
    
    // Write parameters to file
    PrintNumVectorToFile(myfile, parameters);
    
    
    if(i%100==0)
    {
      Rcout << "Iteration " << i << " completed, number added by algorithm = " << Va.length() << std::endl;
    }
    
  }
  
  // Close files
  myfile.close();
  myfile2.close();
  myfile3.close();
  myfile4.close();
  myfile5.close();
  myfile6.close();
  myfile7.close();
  
  double beta_h_prob = (double)nacc_beta_h/(double)iterations;
  double beta_p_prob = (double)nacc_beta_p/(double)iterations;
  double lambda_prob = (double)nacc_lambda/(double)iterations;
  double move_prob = (double)nacc_move/(double)(augmented_moves_proposed[0]);
  double add_prob = (double)nacc_add/(double)(augmented_moves_proposed[1]);
  double remove_prob = (double)nacc_remove/(double)(augmented_moves_proposed[2]);
  double change_prob = (double)nacc_change/(double)(augmented_moves_proposed[3]);
  double swap_prob = (double)nacc_swap/(double)(augmented_moves_proposed[4]);
  double change_random_prob = (double)nacc_change_random/(double)(augmented_moves_proposed[5]);
  double change_all_prob = (double)nacc_change_all/(double)(augmented_moves_proposed[6]);

  Rcout << "beta_P acceptance = " << beta_p_prob << ", beta_H acceptance = " << beta_h_prob << ", lambda acceptance = " << lambda_prob << std::endl;
  Rcout << "Move prob = " << move_prob << ", add prob = " << add_prob << ", remove prob = " << remove_prob 
        << ", change prob = " << change_prob << ", swap prob = " << swap_prob << ", change all prob = " << change_all_prob 
        << ", change random prob = " << change_random_prob << std::endl;
  //Rcout << "Move prob = " << move_prob << ", add prob = " << add_prob << ", remove prob = " << remove_prob << std::endl;
}

