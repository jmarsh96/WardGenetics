#include <Rcpp.h>
#include <iostream>
#include <fstream>
using namespace Rcpp;



IntegerVector WhichVec(int x, IntegerVector vec);
int ReturnSourceSequence_WHD(int sequence_loc, List data);
IntegerVector ReturnDuplicateSequences(int sequence_loc, List data);


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
  if(print_details)
  {
    Rcout << "Screen loglik = " << screening_likelihood << ", import loglik = " << importation_likelihood << ", transmission loglik = " 
          << transmission_likelihood << std::endl;
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
IntegerVector ReturnGenSourceVector_WHD(List data)
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
    int source = ReturnSourceSequence_WHD(i, data);
    gen_source_vector.push_back(source);
  }
  //Rcout << "End calculating gen source\n";
  //Rcout << "After gen source" << std::endl;
  return gen_source_vector;
}

// [[Rcpp::export]]
int ReturnSourceSequence_WHD(int sequence_loc, List data)
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
      int import_sum = WhichVec(-1, data["source"]).length();
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


// [[Rcpp::export]]
void MCMC_EPI_GEN_HCW(List MCMC_options, int sequence_length, IntegerVector t_a, IntegerVector t_c, IntegerVector t_d, IntegerVector source, 
                        IntegerVector hcw_ind, IntegerMatrix screening_matrix, IntegerVector genetic_ids, IntegerVector sample_times,
                        IntegerVector variant_numbers, IntegerMatrix genetic_matrix)
{
  // Initialise data
  List data = InitialiseData_GEN(sequence_length, t_a, t_c, t_d, source, hcw_ind, screening_matrix, genetic_ids, sample_times, 
                                 variant_numbers, genetic_matrix);
  
  
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
      int import_sum = WhichVec(-1, data["source"]).length();
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

