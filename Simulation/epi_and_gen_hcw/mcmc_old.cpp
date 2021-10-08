#include <Rcpp.h>
#include <Rmath.h>
#include <iostream>
#include <fstream>
using namespace Rcpp;
#include <Rcpp/Benchmark/Timer.h>

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::plugins("cpp11")]]



IntegerVector WhichVec(int x, IntegerVector vec);
IntegerVector ReturnGenSourceVector_WHD(List data);
IntegerVector ReturnDuplicatedSequencesTarget_WHD(int target, int time, int current_variant, List data);
IntegerVector ReturnPathToRoot_Rcpp(IntegerVector gen_source, int node);
List UpdateImputedNodes(List data, NumericVector parameters);
bool DoesTargetHaveSwabAtTime(List data, int target, int current_variant, int time);
IntegerVector Rcpp_sort(IntegerVector x, IntegerVector y);
int CalculateDistanceBetweenNodes_Rcpp(int node1, int node2, IntegerVector gen_source, IntegerMatrix genetic_matrix);
IntegerMatrix crossout(IntegerMatrix X, int e);
bool CheckIfThereIsASwabAtTime(List data, int target, int time, int current_variant);
List UpdateImputedNodesMove(List data_cur, List data_can, NumericVector parameters, int target);
List UpdateImputedNodesMove2(List data_cur, List data_can, NumericVector parameters, int target);
List UpdateImputedNodesMove3(List data_cur, List data_can, NumericVector parameters, int target);
void PrintMatrix(IntegerMatrix matrix);
bool DoesNodeNeedToBeImputed(List data, int target, int current_variant);
List UpdateImputedNodesMoveAll(List data_cur, List data_can, NumericVector parameters);
List UniqueSort(IntegerVector x, IntegerVector y, IntegerVector z);
//List UniqueSort(IntegerVector x, IntegerVector y);


List InitialiseData(IntegerVector t_a, IntegerVector t_c, IntegerVector t_d, IntegerVector source, IntegerMatrix screening_matrix)
{
  IntegerVector zero_based_source = clone(source);
  for(int i = 0; i<source.length(); i++)
  {
    if(source[i] > 0)
    {
      zero_based_source[i]--;
    }
  }
  List data = List::create(Named("t_a") = t_a , 
                           Named("t_c") = t_c,
                           Named("t_d") = t_d,
                           Named("source") = zero_based_source,
                           Named("screening_matrix") = screening_matrix);
  return data;
}

// [[Rcpp::export]]
List InitialiseData(int sequence_length, IntegerVector t_a, IntegerVector t_c, IntegerVector t_d, IntegerVector source, IntegerMatrix screening_matrix, 
                    IntegerVector genetic_ids, IntegerVector sample_times, IntegerVector variant_numbers, IntegerMatrix genetic_matrix)
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
                           Named("imputed_nodes") = imputed_nodes);
  return data;
}



void PrintMatrixToFile(std::ofstream &myfile, IntegerMatrix mat)
{
  for(int i = 0; i<mat.nrow(); i++)
  {
    for(int j = 0; j<mat.ncol(); j++)
    {
      myfile << mat(i,j) << " ";
    }
    myfile << std::endl;
  }
  myfile << std::endl;
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

// check if  there is a duplicate sequence at the same time, if so return all sequences which are the duplicate
// [[Rcpp::export]]
IntegerVector ReturnDuplicateSequences2(int sequence_loc, List data)
{
  bool print_debug = false;
  
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
  int inf_time = t_c[patient_id];
  
  
  if(print_debug) 
  {
    Rcout << "Patient ID = " << patient_id << ", sequence time = " << seq_sample_time << ", variant = " << current_variant
          << ", infection time = " << inf_time << "\n";
  }
  
  // also check if there are other swabs at this time (shouldn't happen but will do)
  IntegerVector patient_swab_loc = WhichVec(patient_id, genetic_ids);
  IntegerVector variant_loc = WhichVec(current_variant, variant_numbers);
  IntegerVector patient_variant_swab_locs = intersect(patient_swab_loc, variant_loc);
  IntegerVector patient_swab_times = as<IntegerVector>(sample_times[patient_variant_swab_locs]);
  IntegerVector patient_swab_loc_same_time = WhichVec(seq_sample_time, patient_swab_times);
  if(print_debug)
  {
    Rcout << "Patient swab locs = " << patient_swab_loc << "\n";
    Rcout << "Patient variant swab locs = " << patient_variant_swab_locs << "\n";
    Rcout << "Sequence sample time = " << seq_sample_time << "\n";
    Rcout << "Patient swab times = " << patient_swab_times << "\n";
    Rcout << "Patient swabs location at the same time = " << patient_swab_loc_same_time << "\n";
  }
  
  
  // if there are mu
  if(patient_swab_loc_same_time.length() > 1)
  {
    for(int i = 0; i < patient_swab_loc_same_time.length(); i++)
    {
      int patient_loc = patient_swab_loc_same_time[i];
      duplicated_sequences.push_back(patient_variant_swab_locs[patient_loc]);
    }
  }
  if(print_debug) Rcout << "Duplicated sequences = " << duplicated_sequences << "\n";
  
  if(seq_sample_time == inf_time) 
  {
    // the sequence is observed at the time of infection, therefore we need to consider a source swab at this time, or source infections
    int source = source_vector[patient_id];
    IntegerVector source_swab_loc = WhichVec(source, genetic_ids);
    IntegerVector source_variant_swab_locs = intersect(source_swab_loc, variant_loc);
    
    IntegerVector source_swab_times = as<IntegerVector>(sample_times[source_variant_swab_locs]);
    IntegerVector source_swab_loc_same_time = WhichVec(seq_sample_time, source_swab_times);
    if(print_debug) Rcout << "Source = " << source << ", source swab times = " << source_swab_times << "\n";
    if(source_swab_loc_same_time.length() > 0)
    {
      for(int i = 0; i < source_swab_loc_same_time.length(); i++)
      {
        int source_loc = source_swab_loc_same_time[i];
        duplicated_sequences.push_back(source_variant_swab_locs[source_loc]);
      }
    }
    
    // now check the source infections, if there is a sequence taken at the time of infection, this must also be a duplicate
    IntegerVector source_infections = WhichVec(source, source_vector);
    if(print_debug) Rcout << "Source infections = " << source_infections << "\n";
    for(int i = 0; i < source_infections.length(); i++)
    {
      int target = source_infections[i];
      if(target != patient_id)
      {
        int target_infection_time = t_c[target];
        if(target_infection_time == seq_sample_time)
        {
          // the target was infected at the time of the sequence, so check if they have any swabs at the time
          IntegerVector target_swab_loc = WhichVec(target, genetic_ids);
          IntegerVector target_variant_swab_locs = intersect(target_swab_loc, variant_loc);
          
          IntegerVector target_swab_times = as<IntegerVector>(sample_times[target_variant_swab_locs]);
          IntegerVector target_swab_loc_same_time = WhichVec(seq_sample_time, target_swab_times);
          if(print_debug) Rcout << "Target swabs at the same time = " << target_swab_loc_same_time << "\n";
          if(target_swab_loc_same_time.length() > 0)
          {
            for(int i = 0; i < target_swab_loc_same_time.length(); i++)
            {
              int target_loc = target_swab_loc_same_time[i];
              duplicated_sequences.push_back(target_variant_swab_locs[target_loc]);
              if(print_debug) Rcout << "Duplicated sequences = " << duplicated_sequences << "\n";
            }
          }
        }
      }
    }
  }
  else
  {
    // the sequence is not at the time of infection, therefore we need to look at the current patients infections rather than the source
    IntegerVector source_infections = WhichVec(patient_id, source_vector);
    if(print_debug) Rcout << "Source infections = " << source_infections << "\n";
    for(int i = 0; i < source_infections.length(); i++)
    {
      int target = source_infections[i];
      if(target != patient_id)
      {
        int target_infection_time = t_c[target];
        if(target_infection_time == seq_sample_time)
        {
          // the target was infected at the time of the sequence, so check if they have any swabs at the time
          IntegerVector target_swab_loc = WhichVec(target, genetic_ids);
          IntegerVector target_variant_swab_locs = intersect(target_swab_loc, variant_loc);
          
          IntegerVector target_swab_times = as<IntegerVector>(sample_times[target_variant_swab_locs]);
          IntegerVector target_swab_loc_same_time = WhichVec(seq_sample_time, target_swab_times);
          
          if(target_swab_loc_same_time.length() > 0)
          {
            for(int i = 0; i < target_swab_loc_same_time.length(); i++)
            {
              int target_loc = target_swab_loc_same_time[i];
              duplicated_sequences.push_back(target_variant_swab_locs[target_loc]);
            }
          }
        }
      }
    }
  }
  if(duplicated_sequences.length()>0)
  {
    // check if there are duplicated sequence
    IntegerVector location_in_dup_sequences = WhichVec(sequence_loc, duplicated_sequences);
    if(location_in_dup_sequences.length()==0)
    {
      duplicated_sequences.push_back(sequence_loc);
    }
  }
  return sort_unique(duplicated_sequences);
}


// check if  there is a duplicate sequence at the same time, if so return all sequences which are the duplicate
// NEED TO FIX VARIANTS LATER
// [[Rcpp::export]]
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



int ReturnSourceSequenceTarget_WHD(int target, int time, int current_variant, List data)
{
  IntegerVector genetic_ids = data["genetic_ids"];
  IntegerVector source_vector = data["source"];
  IntegerVector sample_times = data["sample_times"];
  IntegerVector variant_numbers = data["variant_numbers"];
  IntegerVector t_c = data["t_c"];
  
  IntegerVector variant_idx = WhichVec(current_variant, variant_numbers);
  // Check for duplicates at the time
  IntegerVector duplicated_sequences = ReturnDuplicatedSequencesTarget_WHD(target, time, current_variant, data);
  
  
  return 1;
}

// Return a vector of sequence indexs which correspond to swabs or swabs of infectees at the time of infection
IntegerVector ReturnTargetSequences(int target, int current_variant, List data)
{
  IntegerVector genetic_ids = data["genetic_ids"];
  IntegerVector source_vector = data["source"];
  IntegerVector sample_times = data["sample_times"];
  IntegerVector variant_numbers = data["variant_numbers"];
  IntegerVector t_c = data["t_c"];
  
  IntegerVector variant_idx = WhichVec(current_variant, variant_numbers);
  
  IntegerVector patient_swab_loc = WhichVec(target, genetic_ids);
  IntegerVector variant_loc = WhichVec(current_variant, variant_numbers);
  IntegerVector patient_variant_swab_locs = intersect(patient_swab_loc, variant_loc);
  IntegerVector patient_swab_times = as<IntegerVector>(sample_times[patient_variant_swab_locs]);
  //IntegerVector patient_swab_loc_same_time = WhichVec(seq_sample_time, patient_swab_times);
  
  
  // PLACEHOLDER
  return genetic_ids;
}


// Function to return all the sequences at the time in the target, by looking at swabs and infections
// [[Rcpp::export]]
IntegerVector ReturnDuplicatedSequencesTarget_WHD(int target, int time, int current_variant, List data)
{
  IntegerVector genetic_ids = data["genetic_ids"];
  IntegerVector source_vector = data["source"];
  IntegerVector sample_times = data["sample_times"];
  IntegerVector variant_numbers = data["variant_numbers"];
  IntegerVector t_c = data["t_c"];
  
  IntegerVector duplicated_sequences;
  
  // Check if there is any swabs at the time
  IntegerVector patient_swab_loc = WhichVec(target, genetic_ids);
  IntegerVector variant_loc = WhichVec(current_variant, variant_numbers);
  IntegerVector patient_variant_swab_locs = intersect(patient_swab_loc, variant_loc);
  IntegerVector patient_swab_times = as<IntegerVector>(sample_times[patient_variant_swab_locs]);
  IntegerVector patient_swab_loc_same_time = WhichVec(time, patient_swab_times);

  
  
  // Add all sequences at the specified time
  if(patient_swab_loc_same_time.length() > 0)
  {
    for(int i = 0; i < patient_swab_loc_same_time.length(); i++)
    {
      int patient_loc = patient_swab_loc_same_time[i];
      duplicated_sequences.push_back(patient_variant_swab_locs[patient_loc]);
    }
  }
  
  Rcout << "Checked swabs, now infections";
  // Now check all infections
  IntegerVector infections = WhichVec(target, source_vector);
  for(int i = 0; i<infections.length(); i++)
  {
    int current_infectee = infections[i];
    if(t_c[current_infectee] == time)
    {
      // The infectee was infected at the specified time, therefore now check if they have a sequence at that time
      IntegerVector infectee_swab_loc = WhichVec(current_infectee, genetic_ids);
      IntegerVector infectee_variant_swab_locs = intersect(infectee_swab_loc, variant_loc);
      IntegerVector infectee_swab_times = as<IntegerVector>(sample_times[infectee_variant_swab_locs]);
      IntegerVector infectee_swab_loc_same_time = WhichVec(time, infectee_swab_times);
      for(int j = 0; j<patient_swab_loc_same_time.length(); j++)
      {
        int infectee_loc = infectee_swab_loc_same_time[j];
        duplicated_sequences.push_back(infectee_variant_swab_locs[infectee_loc]);
      }
    }
  }
  
  return sort_unique(duplicated_sequences);
}

IntegerMatrix crossout(IntegerMatrix X, int e) {
  IntegerMatrix result = Rcpp::no_init_matrix(X.nrow() - 1, X.ncol() - 1);
  IntegerMatrix::iterator src = X.begin(), dst = result.begin();
  
  while (src != X.end()) {
    if (((src - X.begin()) % X.nrow()) != e && ((src - X.begin()) / X.nrow()) != e) {
      *dst = *src;
      ++dst;
    }
    ++src;
  }
  
  return result;
}

// [[Rcpp::export]]
int ReturnSourceSequence_WHD_old(int sequence_loc, List data)
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
        if(target_inf_time > most_recent_inf_time && target_inf_time < sequence_time)
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
  
  // now a sequence has been found, check if there are duplicates
  IntegerVector source_duplicates = ReturnDuplicateSequences(source_sequence, data);
  if(source_duplicates.length() > 0)
  {
    // there are duplicates of the source sequence, target the last sequence
    if(print_debug) Rcout << "Source duplicates = " << source_duplicates << "\n";
    return source_duplicates[source_duplicates.length()-1];
  }
  
  return source_sequence;
}

int ReturnSourceSequence_WHD_new(int sequence_loc, List data)
{
  bool print_debug = true;
  
  IntegerVector genetic_ids = data["genetic_ids"];
  IntegerVector source_vector = data["source"];
  IntegerVector sample_times = data["sample_times"];
  IntegerVector variant_numbers = data["variant_numbers"];
  IntegerVector t_c = data["t_c"];
  
  // First check if the sequence is a duplicate
  IntegerVector duplicated_sequences = ReturnDuplicateSequences(sequence_loc, data);
  if(duplicated_sequences.length() > 0)
  {
    // There exists duplicate sequences, enforce each one points to the previous one, unless in position zero, then it must look backwards
    if(print_debug) Rcout << "Duplicated sequences = " << duplicated_sequences << "\n";
    int sequence_idx = WhichVec(sequence_loc, duplicated_sequences)[0];
    if(sequence_idx != 0)
    {
      // return the previous one
      if(print_debug)
      {
        Rcout << "Return sequence = " << duplicated_sequences[sequence_idx-1] << " with loc = " << sequence_idx-1 << std::endl;
      }
      return duplicated_sequences[sequence_idx-1];
    }
  }
  
  // now look back to find the most recent swab or infection
  int current_variant = variant_numbers[sequence_loc];
  IntegerVector variant_idx = WhichVec(current_variant, variant_numbers);
  int target = genetic_ids[sequence_loc];
  int sequence_time = sample_times[sequence_loc];
  
  if(sequence_time == t_c[target])
  {
    // The sequence is sampled at the time of infection, therefore we set the target to the source and go from there
    target = source_vector[target];
    if(target==-1)
    {
      // The target is in fact an importation, return -1
      return -1;
    }
  }
  
  //bool sequence_found = false;
  //int source_sequence = -1;
  if(print_debug) Rcout << "Target = " << target << ", sequence time = " << sequence_time << "\n";
  
  /*
   * 
   
  while(!sequence_found)
  {
    
  }
   */
  
  return -1;
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

// [[Rcpp::export]]
int ReturnSourceSequence2_WHD(int sequence_loc, List data)
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
  IntegerVector duplicated_sequences = ReturnDuplicateSequences2(sequence_loc, data);
  if(print_debug) Rcout << "Duplicated sequences = " << duplicated_sequences << "\n";
  if(duplicated_sequences.length() > 0)
  {
    // there exists duplicate sequences, enforce each one points to the previous one, unless in position zero, then it must look backwards
    int sequence_idx = WhichVec(sequence_loc, duplicated_sequences)[0];
    if(sequence_idx != 0)
    {
      // return the previous one
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
  
  // now a sequence has been found, check if there are duplicates
  IntegerVector source_duplicates = ReturnDuplicateSequences2(source_sequence, data);
  if(source_duplicates.length() > 0)
  {
    // there are duplicates of the source sequence, target the last sequence
    if(print_debug) Rcout << "Source duplicates = " << source_duplicates << "\n";
    return source_duplicates[source_duplicates.length()-1];
  }

  return source_sequence;
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
IntegerVector ReturnGenSourceVector2_WHD(List data)
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
    int source = ReturnSourceSequence2_WHD(i, data);
    gen_source_vector.push_back(source);
  }
  //Rcout << "End calculating gen source\n";
  //Rcout << "After gen source" << std::endl;
  return gen_source_vector;
}


// this function will take the current configuration of the tree and return a new tree /genetic ids that need to be imputed
// to capture the shared stochastic evolution of the sequences
// for now assume that there are no nodes that are imputed, however later the function will have to compute the difference
// [[Rcpp::export]]
List ReturnUpdatedData(List input_data)
{
  List data = clone(input_data);
  IntegerVector genetic_ids = data["genetic_ids"];
  IntegerVector sample_times = data["sample_times"];
  IntegerVector t_c = data["t_c"];
  IntegerVector source = data["source"];
  IntegerVector variant_numbers = data["variant_numbers"];
  
  int original_length = genetic_ids.length();
  IntegerVector updated_genetic_ids = clone(genetic_ids);
  IntegerVector updated_sample_times = clone(sample_times);
  
  for(int i = 0; i<t_c.length(); i++)
  {
    if(t_c[i] != -1 && source[i] != -1)
    {
      // check if there is a sequence at the time of infection
      int infection_time = t_c[i];
      IntegerVector gen_locs = WhichVec(i, genetic_ids);
      IntegerVector gen_times = as<IntegerVector>(sample_times[gen_locs]);
      IntegerVector inf_loc = WhichVec(infection_time, gen_times);
      if(inf_loc.length() == 0)
      {
        // there is not a sequence at the time of infection, therefore add one
        genetic_ids.push_back(i);
        sample_times.push_back(t_c[i]);
      }
    }
  }
  
  Rcout << "After adding nodes, genetic ids = " << genetic_ids << std::endl;
  Rcout << "After adding nodes, sample times = " << sample_times << std::endl;
  
  List new_data = List::create(Named("t_c") = t_c,
                           Named("source") = source,
                           Named("genetic_ids") = genetic_ids,
                           Named("sample_times") = sample_times,
                           Named("variant_numbers") = variant_numbers);
  
  IntegerVector updated_gen_source_vector = ReturnGenSourceVector_WHD(new_data);
  Rcout << "Updated gen source vector = " << updated_gen_source_vector << std::endl;
  for(int i = original_length; i<updated_gen_source_vector.length(); i++)
  {
    int current_id = genetic_ids[i];
    IntegerVector gen_source_locs = WhichVec(i, updated_gen_source_vector);
    Rcout << "Current id = " << current_id << ", gen source locs = " << gen_source_locs << std::endl;
    if(gen_source_locs.length() > 1)
    {
      // keep the node, it is helpful
      updated_genetic_ids.push_back(genetic_ids[i]);
      updated_sample_times.push_back(sample_times[i]);
    }
  }
  
  Rcout << "Updated genetic ids = " << updated_genetic_ids << std::endl;
  Rcout << "Updated sample times = " << updated_sample_times << std::endl;
  
  return data;
}


// [[Rcpp::export]]
int col_pop(List data, int day)
{
  int cp = 0;
  IntegerVector t_c = data["t_c"];
  IntegerVector t_d = data["t_d"];
  IntegerVector source = data["source"];
  for(int i = 0; i<t_c.length(); i++)
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




int CalculateTruePositives(List data)
{
  IntegerMatrix screening_matrix = data["screening_matrix"];
  int TP = 0;
  for(int i = 0; i<screening_matrix.nrow(); i++)
  {
    IntegerVector y = WhichVec(1,screening_matrix(i,_));
    TP += y.length();
  }
  return TP;
}

/*

// [[Rcpp::export]]
List UniqueSort(IntegerVector x, IntegerVector y)
{
  // vector x must already be sorted by y in ascending order
  

  int i = 0;
  IntegerVector x_out(x.length());
  while(i < y.length())
  {
    int current_time = y[i];
    
    // check if there are any other times the same
    IntegerVector same_time_idx = WhichVec(current_time, y);
    if(same_time_idx.length()==1)
    {
      x_out[i] = x[i];
      i++;
    }
    else
    {
      // Other times have been found, need to order them uniquely
      IntegerVector x_ids = as<IntegerVector>(x[same_time_idx]);
      IntegerVector new_time_idx = Rcpp_sort(same_time_idx, x_ids);
      for(int j = 0; j<new_time_idx.length(); j++)
      {
        int old_idx = i + j;
        int new_idx = new_time_idx[j];
        x_out[new_idx] = x[old_idx];
        //Rcout << "old idx = " << old_idx << ", new idx = " << new_idx << std::endl;
      }
      i += x_ids.length();
    }
  }
  List final_data = List::create(Named("nodes") = x_out, 
                                 Named("times") = y);
  return final_data;
}

*/
 
// [[Rcpp::export]]
List UniqueSort(IntegerVector x, IntegerVector y, IntegerVector z)
{
  // vector x must already be sorted by y in ascending order
  
  
  int i = 0;
  IntegerVector x_out(x.length());
  IntegerVector z_out(z.length());
  while(i < y.length())
  {
    int current_time = y[i];
    
    // check if there are any other times the same
    IntegerVector same_time_idx = WhichVec(current_time, y);
    if(same_time_idx.length()==1)
    {
      x_out[i] = x[i];
      z_out[i] = z[i];
      i++;
    }
    else
    {
      // Other times have been found, need to order them uniquely
      IntegerVector x_ids = as<IntegerVector>(x[same_time_idx]);
      IntegerVector new_time_idx = Rcpp_sort(same_time_idx, x_ids);
      for(int j = 0; j<new_time_idx.length(); j++)
      {
        int old_idx = i + j;
        int new_idx = new_time_idx[j];
        x_out[new_idx] = x[old_idx];
        z_out[new_idx] = z[old_idx];
        //Rcout << "old idx = " << old_idx << ", new idx = " << new_idx << std::endl;
      }
      i += x_ids.length();
    }
  }
  List final_data = List::create(Named("nodes") = x_out, 
                                 Named("times") = y,
                                 Named("variant_numbers") = z_out);
  return final_data;
}

// [[Rcpp::export]]
int CalculateFalseNegatives(List data)
{
  IntegerVector t_c = data["t_c"];
  IntegerVector t_d = data["t_d"];
  IntegerMatrix screening_matrix = data["screening_matrix"];
  int FN = 0;
  for(int i = 0; i<t_c.length(); i++)
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

int WhichMin(IntegerVector vec)
{
  int length = vec.length();
  if(length == 0) stop("empty vector to find minimum element");
  if(length==1)
  {
    return 0;
  }
  else
  {
    int min = vec[0];
    int min_loc = 0;
    for(int i = 1; i<length; i++)
    {
      int current = vec[i];
      if(current < min)
      {
        min = current;
        min_loc = i;
      }
    }
    return min_loc;
  }
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



double CalculateTransmissionLikelihood(List data, NumericVector parameters)
{
  double loglik = 0.0;
  IntegerVector t_a = data["t_a"];
  IntegerVector t_c = data["t_c"];
  IntegerVector t_d = data["t_d"];
  IntegerVector source = data["source"];
  double beta = parameters[2];
  
  for(int i = 0; i<t_c.length(); i++)
  {
    if(t_c[i] == -1)
    {
      // The individual avoided infection the whole time
      for(int t = t_a[i]; t<=t_d[i]; t++)
      {
        loglik -= beta*col_pop(data,t);
      }
    }
    else if(source[i] == -1)
    {
      // The individual is an importation, and so no contribution to the likelihood
    }
    else if(t_c[i] == t_a[i])
    {
      // The individual was infected on their first day
      loglik += log(1-exp(-beta*col_pop(data,t_a[i]))) - log(col_pop(data,t_c[i]));
    }
    else
    {
      // The individual was infected after the first day
      for(int t = t_a[i]; t<t_c[i]; t++)
      {
        loglik -= beta*col_pop(data,t);
      }
      loglik += log(1-exp(-beta*col_pop(data,t_c[i]))) - log(col_pop(data,t_c[i]));
    }
    //Rcout << "i = " << i << ", loglik = " << loglik << ", col pop = " << col_pop(data,t_c[i]) << ", source = " << source[i] << std::endl;
    
    // perform a quick check to make sure their infector is actually able to infect them
    if(t_c[i] != -1)
    {
      int infector = source[i];
      if(infector != -1)
      {
        // ignore if they are an importation
        int infector_source = source[infector];
        if(infector_source == -1)
        {
          // infector is an import
          if(!(t_c[infector] <= t_c[i] && t_d[infector] >= t_c[i] && t_c[infector] != -1))
          {
            Rcout << "Target = " << i << ", t_c = " << t_c[i] << ", (import) infector = " 
                  << infector << ", infector t_c = " << t_c[infector] << std::endl;
            stop("Error with likelihood calculation - there is a break in the chain and a source is unable to infect");
          }
        }
        else
        {
          // infector is an acquisition
          if(!(t_c[infector] < t_c[i] && t_d[infector] >= t_c[i] && t_c[infector] != -1))
          {
            Rcout << "Target = " << i << ", t_c = " << t_c[i] << ", infector = " 
                  << infector << ", infector t_c = " << t_c[infector] << std::endl;
            stop("Error with likelihood calculation - there is a break in the chain and a source is unable to infect");
          }
        }
      }

    }
    
  }
  return loglik;
}

double CalculateTransmissionLikelihood_NS(List data, NumericVector parameters)
{
  double loglik = 0.0;
  IntegerVector t_a = data["t_a"];
  IntegerVector t_c = data["t_c"];
  IntegerVector t_d = data["t_d"];
  IntegerVector source = data["source"];
  double beta = parameters[2];
  
  for(int i = 0; i<t_c.length(); i++)
  {
    if(t_c[i] == -1)
    {
      // The individual avoided infection the whole time
      for(int t = t_a[i]; t<=t_d[i]; t++)
      {
        loglik -= beta*col_pop(data,t);
      }
    }
    else if(source[i] == -1)
    {
      // The individual is an importation, and so no contribution to the likelihood
    }
    else if(t_c[i] == t_a[i])
    {
      // The individual was infected on their first day
      loglik += log(1-exp(-beta*col_pop(data,t_a[i])));
    }
    else
    {
      // The individual was infected after the first day
      for(int t = t_a[i]; t<t_c[i]; t++)
      {
        loglik -= beta*col_pop(data,t);
      }
      loglik += log(1-exp(-beta*col_pop(data,t_c[i])));
    }
    //Rcout << "i = " << i << ", loglik = " << loglik << ", col pop = " << col_pop(data,t_c[i]) << std::endl;
    
  }
  return loglik;
}

double CalculateScreeningLikelihood(List data, NumericVector parameters)
{
  double z = parameters[0];
  double loglik = CalculateTruePositives(data) * log(z) + CalculateFalseNegatives(data) * log(1-z);
  return loglik;
}

double CalculateImportationLikelihood(List data, NumericVector parameters)
{
  double p = parameters[1];
  IntegerVector source = data["source"];
  IntegerVector import_idx = WhichVec(-1,source);
  int population_size = source.length();
  int importation_sum = import_idx.length();
  double loglik = importation_sum * log(p) + (population_size-importation_sum) * log(1-p);
  return loglik;
}

double CalculateMutationProbability(double rate, double t)
{
  double prob = 0.75*(1-exp(-4*rate*t));
  return prob;
}



double CalculateGeneticLikelihood(List data, NumericVector parameters)
{
  double lambda = parameters[3];
  double importation_distance = parameters[4];
  //double variant_distance = parameters[5];
  IntegerVector gen_source = ReturnGenSourceVector_WHD(data);
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

  //Rcout << "Zero time counter = " << zero_time_counter << ", import length = " << import_idx.length() << std::endl;
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
  if(data.containsElementNamed("genetic_ids"))
  {
    loglik += CalculateGeneticLikelihood(data, parameters);
  }
  if(print_details)
  {
    Rcout << "Screen loglik = " << screening_likelihood << ", import loglik = " << importation_likelihood << ", transmission loglik = " 
          << transmission_likelihood << ", genetic loglik = " << CalculateGeneticLikelihood(data, parameters) << std::endl;
  }
  

  
  return loglik;
}

// [[Rcpp::export]]
double CalculateLogLikelihood_NS(List data, NumericVector parameters)
{
  double screening_likelihood = CalculateScreeningLikelihood(data, parameters);
  double importation_likelihood = CalculateImportationLikelihood(data, parameters);
  double transmission_likelihood = CalculateTransmissionLikelihood_NS(data, parameters);
  //Rcout << "Screen loglik = " << screening_likelihood << ", import loglik = " << importation_likelihood << ", transmission loglik = " << transmission_likelihood << std::endl;
  double loglik = screening_likelihood + importation_likelihood + transmission_likelihood;
  return loglik;
}

int CalculateImportSum(List data)
{
  IntegerVector source = data["source"];
  IntegerVector import_locs = WhichVec(-1,source);
  return import_locs.length();
}

void UpdateTransmissionRate(List data, NumericVector &parameters, double &loglik_cur, double proposal_variance, double prior, int &nacc)
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

void UpdateMutationRate(List data, NumericVector &parameters, double &loglik_cur, double proposal_variance, double prior, int &nacc)
{
  NumericVector parameters_can = clone(parameters);
  parameters_can[3] = R::rnorm(parameters[3], proposal_variance);
  if(parameters_can[3] > 0)
  {
    double logpi_cur = loglik_cur + R::dexp(parameters[3], 1/prior, 1);
    double loglik_can = CalculateLogLikelihood(data, parameters_can);
    double logpi_can = loglik_can + R::dexp(parameters_can[3], 1/prior, 1);
    double U = R::runif(0.0,1.0);
    if(log(U) < logpi_can - logpi_cur)
    {
      parameters = parameters_can;
      loglik_cur = loglik_can;
      nacc++;
    }
  }
}

void UpdateTransmissionRate_NS(List data, NumericVector &parameters, double &loglik_cur, double proposal_variance, double prior, int &nacc)
{
  NumericVector parameters_can = clone(parameters);
  parameters_can[2] = R::rnorm(parameters[2], proposal_variance);
  if(parameters_can[2] > 0)
  {
    double logpi_cur = loglik_cur + R::dexp(parameters[2], 1/prior, 1);
    double loglik_can = CalculateLogLikelihood_NS(data, parameters_can);
    double logpi_can = loglik_can + R::dexp(parameters_can[2], 1/prior, 1);
    double U = R::runif(0.0,1.0);
    if(log(U) < logpi_can - logpi_cur)
    {
      parameters = parameters_can;
      loglik_cur = loglik_can;
      nacc++;
    }
  }
}

// Vq is the set of patients who tested positive during their stay
IntegerVector CalculateVq(List data)
{
  IntegerVector Vq;
  IntegerMatrix screening_matrix = data["screening_matrix"];
  for(int i = 0; i<screening_matrix.nrow(); i++)
  {
    IntegerVector current_row = screening_matrix(i,_);
    IntegerVector positive_locs = WhichVec(1, current_row);
    if(positive_locs.length()>0)
    {
      Vq.push_back(i);
    }
  }
  return Vq;
}

IntegerVector CalculateVs(List data)
{
  IntegerVector Vs;
  IntegerMatrix screening_matrix = data["screening_matrix"];
  for(int i = 0; i<screening_matrix.nrow(); i++)
  {
    IntegerVector current_row = screening_matrix(i,_);
    IntegerVector positive_locs = WhichVec(1, current_row);
    if(positive_locs.length()==0)
    {
      Vs.push_back(i);
    }
  }
  return Vs;
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



int SampleVector(IntegerVector x)
{
  double U = R::runif(0.0,1.0);
  int length = x.length();
  return x[floor(U*length)];
}



int ReturnLastDay(List data, int target, int import_ind)
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
    int offspring_coltime = t_c[offspring[i]] - 1;
    //Rcout << offspring_coltime << " ";
    if(offspring_coltime < last_day)
    {
      if(import_ind == 1)
      {
        last_day = offspring_coltime; // not minus one beacuse they are an import
      }
      else
      {
        last_day = offspring_coltime;// - 1;  // minus one because you can only infect someone on the day before
      }
    }
  }
  
  IntegerVector screening_row = screening_matrix(target,_);
  IntegerVector positive_days = WhichVec(1,screening_row);
  //Rcout << "positive days = " << positive_days;
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

IntegerVector ReturnPossibleInfectors(List data, int target)
{
  IntegerVector t_c = data["t_c"];
  IntegerVector t_d = data["t_d"];
  IntegerVector source = data["source"];
  IntegerVector possible_infectors;
  
  int day_of_colonisation = t_c[target];
  for(int i = 0; i<t_c.length(); i++)
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



// check everything is ok
bool CheckData(List data)
{
  IntegerVector t_c = data["t_c"];
  IntegerVector t_d = data["t_d"];
  IntegerVector source = data["source"];
  for(int i = 0; i<t_c.length(); i++)
  {
    if(t_c[i] != -1)
    {
      int current_source = source[i];
      if(current_source != -1)
      {
        int sources_source = source[current_source];
        if(sources_source == -1)
        {
          // The current source is an importation, therefore check they are in the ward and able to infect
          if(t_c[i] >= t_c[current_source] && t_c[i] <= t_d[current_source])
          {
            
          }
          else
          {
            Rcout << "Error with ID = " << i << ", source = " << current_source << std::endl;
            return false;
          }
        }
        else
        {
          // the current source is an acquisiton, therefore check there are in ward and able to infect
          if(t_c[i] > t_c[current_source] && t_c[i] <= t_d[current_source])
          {
            
          }
          else
          {
            return false;
            Rcout << "Error with ID = " << i << ", source = " << current_source << std::endl;
          }
        }
      }
    }
  }
  return true;
}


double CalculateGeneticContributionMove(List data_cur, List &data_can)
{
  IntegerVector gen_source_cur = ReturnGenSourceVector_WHD(data_cur);
  IntegerVector gen_source_can = ReturnGenSourceVector_WHD(data_can);
  IntegerVector imputed_distances_can = data_can["imputed_distances"];
  int sequence_length = data_can["sequence_length"];
  double importation_p = data_can["importation_p"];
  double contribution = 0;
  for(int i = 0; i<gen_source_cur.length(); i++)
  {
    if(gen_source_cur[i] == -1 && gen_source_can[i] != -1)
    {
      // a previous import is no longer, so remove the imputed distance
      contribution += R::dbinom(imputed_distances_can[i], sequence_length, importation_p, 1);
      imputed_distances_can[i] = 0;
    }
    else if(gen_source_cur[i] != -1 && gen_source_can[i] == -1)
    {
      // the proposed move has an imported idx, so propose a new distance
      int draw = R::rbinom(sequence_length, importation_p);
      imputed_distances_can[i] = draw;
      contribution -= R::dbinom(imputed_distances_can[i], sequence_length, importation_p, 1);
    }
  }
  return contribution;
}




void MoveColonisationTime(List &data, NumericVector parameters, double &loglik_cur, double w, IntegerVector Vq, IntegerVector Va, int &nacc_move)
{
  bool print_output = false;
  IntegerVector t_a = data["t_a"];
  IntegerVector t_c = data["t_c"];
  IntegerVector t_d = data["t_d"];
  IntegerVector source = data["source"];
  
  List data_can = clone(data);
  IntegerVector t_c_can = data_can["t_c"];
  IntegerVector source_can = data_can["source"];
  
  IntegerVector colonised_population = union_(Vq,Va);
  int target = SampleVector(colonised_population);
  
  //Rcout << "target = " << target << ", last day = " << last_day << std::endl;
  double log_prop_ratio = 0;
  
  double U = R::runif(0.0,1.0);
  if(U < w)
  {
    // Propose the target is an importation
    t_c_can[target] = t_a[target];
    source_can[target] = -1;
    if(source[target] == -1)
    {
      // Importation -> Importation
      // Do nothing, proposal ratio = 1
      
      if(print_output) Rcout << "Move Imp -> Imp";
    }
    else
    {
      // Acquisition -> Importation
      int last_day = ReturnLastDay(data, target, 0);
      log_prop_ratio = log(1-w) - log(w) - log(last_day-t_a[target]+1) - log(col_pop(data, t_c[target]));
      if(print_output) Rcout << "Move Acq -> Imp";
    }
  }
  else
  {
    // Propose the target is an acquisition
    //Rcout << "t_a = " << t_a[target] << ", last day = " << last_day << std::endl;
    int last_day = ReturnLastDay(data, target, 0);
    if(last_day == (t_a[target] - 1))
    {
      // the person is an importation and has infected someone on day one, therefore it is not possible
      // to reassign them as an acquisition
      return;
    }
    t_c_can[target] = SampleVector(seq(t_a[target],last_day));
    IntegerVector possible_infectors = ReturnPossibleInfectors(data_can, target);
    if(possible_infectors.length()==0) return;
    source_can[target] = SampleVector(possible_infectors);
    if(possible_infectors.length() != col_pop(data_can, t_c_can[target]))
    {
      Rcout << "Target = " << target << ", Possible infectors can = " << possible_infectors << std::endl;
      Rcout << "C(t*) = " << col_pop(data_can, t_c_can[target]) << std::endl;
      Rcout << "t_a = " << t_a[target] << ", t_d = " << t_d[target] << ", last_day = " << last_day 
            << ", proposed coltime = " << t_c_can[target] << ", source before = " << source[target] <<", source = " << source_can[target] << std::endl;
      stop("error, candidate infectors are the wrong length");
    }

    if(source[target] != -1)
    {
      if(print_output) Rcout << "Move Acq -> Acq";
      // Acquisition -> Acquisition
      log_prop_ratio = log(col_pop(data_can, t_c_can[target])) - log(col_pop(data, t_c[target]));
      //Rcout << "C(t*) = " << col_pop(data_can, t_c_can[target]) << ", C(t) = " << col_pop(data, t_c[target]) << std::endl;
      IntegerVector possible_infectors_cur = ReturnPossibleInfectors(data, target);
      if(possible_infectors_cur.length() != col_pop(data, t_c[target]))
      {
        stop("error, current infectors are the wrong length");
      }
    }
    else
    {
      if(print_output) Rcout << "Move Imp -> Acq";
      // Importation -> Acquisition
      //log_prop_ratio = log(w) + log(last_day - t_a[target] + 1) + log(col_pop(data_can,t_c_can[target])) - log(1-w);
      log_prop_ratio = log(w) + log(last_day - t_a[target] + 1) + log(col_pop(data_can,t_c_can[target])) - log(1-w);
      
      // check there is no one at the time of infection
    }
  }
  
  //if(!CheckData(data_can)) stop("inconsistency in data");
  bool impute_distances = false;
  if(data.containsElementNamed("genetic_ids") && impute_distances)
  {
    data = UpdateImputedNodes(data, parameters);
    double genetic_contribution = data["genetic_contribution"];
    log_prop_ratio += genetic_contribution;
  }
  //data = UpdateImputedNodes(data, parameters);

  
  double loglik_can = CalculateLogLikelihood(data_can, parameters);
  U = R::runif(0.0,1.0);
  if(print_output) Rcout << ", target = " << target <<  ": Loglik can = " << loglik_can << ", loglik cur = " << loglik_cur << ", log prop ratio = " << log_prop_ratio;
  //Rcout << "t_a = " << t_a[target] << ", t_d = " << t_d[target] << ", last_day = " << last_day 
  //      << ", proposed coltime = " << t_c_can[target] << ", source before = " << source[target] <<", source = " << source_can[target] << std::endl;
  //Rcout << "  MOVE ID " << target <<  ": Loglik can = " << loglik_can << ", loglik cur = " << loglik_cur << ", log prop ratio = " << log_prop_ratio;
  //Rcout << "parameters = " << parameters << std::endl;
  if(log(U) < loglik_can - loglik_cur + log_prop_ratio)
  {
    data = data_can;
    loglik_cur = loglik_can;
    if(print_output) Rcout << " - accepted." << std::endl;
    //Rcout << "Gen source cur = " << ReturnGenSourceVector_WHD(data) << std::endl;
    //Rcout << "Gen source can = " << ReturnGenSourceVector_WHD(data_can) << std::endl;
    nacc_move++;
  }
  else
  {
    if(print_output) Rcout << " - rejected." << std::endl;
    //Rcout << "Gen source cur = " << ReturnGenSourceVector_WHD(data) << std::endl;
    //Rcout << "Gen source can = " << ReturnGenSourceVector_WHD(data_can) << std::endl;
  }
}


List MoveColonisationTimeGenetic(List data, NumericVector parameters, double &loglik_cur, double w, IntegerVector Vq, IntegerVector Va, int &nacc_move)
{
  bool print_output = false;
  bool print_output2 = false;
  IntegerVector t_a = data["t_a"];
  IntegerVector t_c = data["t_c"];
  IntegerVector t_d = data["t_d"];
  IntegerVector source = data["source"];
  
  List data_can = clone(data);
  IntegerVector t_c_can = data_can["t_c"];
  IntegerVector source_can = data_can["source"];
  
  IntegerVector colonised_population = union_(Vq,Va);
  int target = SampleVector(colonised_population);
  
  //Rcout << "target = " << target << ", last day = " << last_day << std::endl;
  double log_prop_ratio = 0;
  if(print_output) Rcout << "target = " << target << std::endl;
  double U = R::runif(0.0,1.0);
  if(U < w)
  {
    // Propose the target is an importation
    t_c_can[target] = t_a[target];
    source_can[target] = -1;
    if(source[target] == -1)
    {
      // Importation -> Importation
      // Do nothing, proposal ratio = 1
      
      if(print_output) Rcout << "Move Imp -> Imp";
    }
    else
    {
      // Acquisition -> Importation
      int last_day = ReturnLastDay(data, target, 0);
      log_prop_ratio = log(1-w) - log(w) - log(last_day-t_a[target]+1) - log(col_pop(data, t_c[target]));
      if(print_output) Rcout << "Move Acq -> Imp";
    }
  }
  else
  {
    // Propose the target is an acquisition
    //Rcout << "t_a = " << t_a[target] << ", last day = " << last_day << std::endl;
    int last_day = ReturnLastDay(data, target, 0);
    if(print_output) Rcout << "t_a = " << t_a[target] << ", last day = " << last_day << std::endl;
    if(last_day == (t_a[target] - 1))
    {
      // the person is an importation and has infected someone on day one, therefore it is not possible
      // to reassign them as an acquisition
      return data;
    }
    t_c_can[target] = SampleVector(seq(t_a[target],last_day));
    IntegerVector possible_infectors = ReturnPossibleInfectors(data_can, target);
    if(print_output) Rcout << "Possible infectors = " << possible_infectors << std::endl;
    if(possible_infectors.length()==0) return data;
    source_can[target] = SampleVector(possible_infectors);
    if(possible_infectors.length() != col_pop(data_can, t_c_can[target]))
    {
      Rcout << "Target = " << target << ", Possible infectors can = " << possible_infectors << std::endl;
      Rcout << "C(t*) = " << col_pop(data_can, t_c_can[target]) << std::endl;
      Rcout << "t_a = " << t_a[target] << ", t_d = " << t_d[target] << ", last_day = " << last_day 
            << ", proposed coltime = " << t_c_can[target] << ", source before = " << source[target] <<", source = " << source_can[target] << std::endl;
      stop("error, candidate infectors are the wrong length");
    }
    
    if(source[target] != -1)
    {
      if(print_output) Rcout << "Move Acq -> Acq";
      // Acquisition -> Acquisition
      log_prop_ratio = log(col_pop(data_can, t_c_can[target])) - log(col_pop(data, t_c[target]));
      //Rcout << "C(t*) = " << col_pop(data_can, t_c_can[target]) << ", C(t) = " << col_pop(data, t_c[target]) << std::endl;
      IntegerVector possible_infectors_cur = ReturnPossibleInfectors(data, target);
      if(possible_infectors_cur.length() != col_pop(data, t_c[target]))
      {
        stop("error, current infectors are the wrong length");
      }
    }
    else
    {
      if(print_output) Rcout << "Move Imp -> Acq";
      // Importation -> Acquisition
      //log_prop_ratio = log(w) + log(last_day - t_a[target] + 1) + log(col_pop(data_can,t_c_can[target])) - log(1-w);
      log_prop_ratio = log(w) + log(last_day - t_a[target] + 1) + log(col_pop(data_can,t_c_can[target])) - log(1-w);
      
      // check there is no one at the time of infection
    }
  }
  
  if(print_output2)
  {
    Rcout << "Propose to move target " << target << ", from t=" << t_c[target] << " to t=" << t_c_can[target] << ", from source " << source[target] << " to "
          << source_can[target] << std::endl;
  }
  
  //if(!CheckData(data_can)) stop("inconsistency in data");
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

  //data = UpdateImputedNodes(data, parameters);
  //Rcout << "Updated imputed nodes " << std::endl;
  
  


  
  double loglik_can = CalculateLogLikelihood(data_can, parameters);
  
  bool specific_check = false;
  if(specific_check && loglik_can < -1000)
  {
    //Rcout << std::endl << "t_c = " << t_c << std::endl;
    //Rcout << "t_c_can = " << t_c_can << std::endl;
    Rcout << std::endl;
    IntegerVector gen_source = ReturnGenSourceVector_WHD(data_can);
    IntegerVector gen_source_cur = ReturnGenSourceVector_WHD(data);
    Rcout << "gen source cur = " << gen_source_cur << std::endl;
    Rcout << "gen source can = " << gen_source << std::endl;
    //stop("return early, problem with likelihood!");
  }
  
  if(print_output2)
  {
    IntegerVector genetic_ids_can = data_can["genetic_ids"];
    IntegerVector sample_times_can = data_can["sample_times"];
    IntegerMatrix genetic_matrix_can = data_can["genetic_matrix"];
    Rcout << "Genetic ids  = " << genetic_ids_can << std::endl;
    Rcout << "Sample times = " << sample_times_can << std::endl;
    Rcout << "Gen source = " << ReturnGenSourceVector_WHD(data_can) << std::endl;
    Rcout << "--- Genetic matrix ---" << std::endl;
    PrintMatrix(genetic_matrix_can);
    Rcout << "------------------" << std::endl;
    Rcout << "ll can = " << loglik_can << ", ll cur = " << loglik_cur << ", log prop ratio = " << log_prop_ratio << std::endl;
    double screening_likelihood = CalculateScreeningLikelihood(data, parameters);
    double importation_likelihood = CalculateImportationLikelihood(data, parameters);
    double transmission_likelihood = CalculateTransmissionLikelihood(data, parameters);
    double genetic_likelihood = CalculateGeneticLikelihood(data, parameters);
    Rcout << "screen ll = " << screening_likelihood << ", import ll = " << importation_likelihood << ", transmission ll = " << transmission_likelihood << ", gen ll = "
          << genetic_likelihood << std::endl;
    double screening_likelihood_can = CalculateScreeningLikelihood(data_can, parameters);
    double importation_likelihood_can = CalculateImportationLikelihood(data_can, parameters);
    double transmission_likelihood_can = CalculateTransmissionLikelihood(data_can, parameters);
    double genetic_likelihood_can = CalculateGeneticLikelihood(data_can, parameters);
    Rcout << "screen ll can = " << screening_likelihood_can << ", import ll can = " << importation_likelihood_can << ", transmission ll can = " << transmission_likelihood_can
          << ", gen ll can = " << genetic_likelihood_can << std::endl;
  }
  //stop("stop early !");
  
  if(!(log_prop_ratio < 1000000 && log_prop_ratio > -1000000))
  {
    Rcout << "LPR = " << log_prop_ratio << std::endl;
    stop("LPR wrong");
  }
  
  U = R::runif(0.0,1.0);
  if(print_output) Rcout << ", target = " << target <<  ": Loglik can = " << loglik_can << ", loglik cur = " << loglik_cur << ", log prop ratio = " << log_prop_ratio;
  //Rcout << "t_a = " << t_a[target] << ", t_d = " << t_d[target] << ", last_day = " << last_day 
  //      << ", proposed coltime = " << t_c_can[target] << ", source before = " << source[target] <<", source = " << source_can[target] << std::endl;
  //Rcout << "  MOVE ID " << target <<  ": Loglik can = " << loglik_can << ", loglik cur = " << loglik_cur << ", log prop ratio = " << log_prop_ratio;
  //Rcout << "parameters = " << parameters << std::endl;
  if(log(U) < loglik_can - loglik_cur + log_prop_ratio)
  {
    //data = data_can;
    loglik_cur = loglik_can;
    if(print_output2) Rcout << "***** Accepted. *****" << std::endl << std::endl;
    if(print_output)
    {
      Rcout << " - accepted." << std::endl;
      IntegerVector gen_source = ReturnGenSourceVector_WHD(data);
      Rcout << "gen source cur = " << gen_source << std::endl;
    }
    //Rcout << "Gen source cur = " << ReturnGenSourceVector_WHD(data) << std::endl;
    //Rcout << "Gen source can = " << ReturnGenSourceVector_WHD(data_can) << std::endl;
    nacc_move++;
    return data_can;
  }
  else
  {
    if(print_output2) Rcout << "***** Rejected. *****" << std::endl << std::endl;
    if(print_output)
    {
      Rcout << " - rejected." << std::endl;
      IntegerVector gen_source = ReturnGenSourceVector_WHD(data);
      Rcout << "gen source cur = " << gen_source << std::endl;
    }
    //Rcout << "Gen source cur = " << ReturnGenSourceVector_WHD(data) << std::endl;
    //Rcout << "Gen source can = " << ReturnGenSourceVector_WHD(data_can) << std::endl;
    return data;
  }
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
    //t_c_can[target] = floor(R::runif(t_a[target],t_d[target]+1));
    IntegerVector possible_infectors = ReturnPossibleInfectors(data_can, target);
    if(possible_infectors.length()==0) return;
    source_can[target] = SampleVector(possible_infectors);
    
    log_prop_ratio = log(col_pop(data_can, t_c_can[target])) + log(susceptible_population.length()) + log(t_d[target] - t_a[target] + 1) - log(1-w) - log(1+nadd0.length());
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

List AddColonisationTime_List(List data, NumericVector parameters, double &loglik_cur, double w, IntegerVector Vs, IntegerVector &Va, int &nacc_add)
{
  IntegerVector t_a = data["t_a"];
  IntegerVector t_c = data["t_c"];
  IntegerVector t_d = data["t_d"];
  IntegerVector source = data["source"];
  
  List data_can = clone(data);
  IntegerVector t_c_can = data_can["t_c"];
  IntegerVector source_can = data_can["source"];
  
  IntegerVector susceptible_population = setdiff(Vs,Va);
  if(susceptible_population.length()==0) return data;
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
    //t_c_can[target] = floor(R::runif(t_a[target],t_d[target]+1));
    IntegerVector possible_infectors = ReturnPossibleInfectors(data_can, target);
    if(possible_infectors.length()==0) return data;
    source_can[target] = SampleVector(possible_infectors);
    
    log_prop_ratio = log(col_pop(data_can, t_c_can[target])) + log(susceptible_population.length()) + log(t_d[target] - t_a[target] + 1) - log(1-w) - log(1+nadd0.length());
  }
  
  bool impute_distances = false;
  if(data.containsElementNamed("genetic_ids") && impute_distances)
  {
    Rcout << "Before update imputed nodes" << std::endl;
    data_can = UpdateImputedNodes(data_can, parameters);
    Rcout << "After update imputed nodes" << std::endl;
    double genetic_contribution = data_can["genetic_contribution"];
    //if(print_output) Rcout << "Gen contribution = " << genetic_contribution << std::endl;
    log_prop_ratio += genetic_contribution;
  }
  
  //if(!CheckData(data_can)) stop("error with data with the adding step!");
  double loglik_can = CalculateLogLikelihood(data_can, parameters);
  U = R::runif(0.0,1.0);
  //Rcout << "ADDING ID " << target <<  ": Loglik can = " << loglik_can << ", loglik cur = " << loglik_cur << ", log prop ratio = " << log_prop_ratio;
  if(log(U) < loglik_can - loglik_cur + log_prop_ratio)
  {
    loglik_cur = loglik_can;
    Va.push_back(target);
    nacc_add++;
    return data_can;
    //Rcout << " - accepted. Number added = " << Va.length() << std::endl;
  }
  else
  {
    return data;
    //Rcout << " - rejected." << std::endl;
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
    log_prop_ratio = log(nadd0.length()) + log(1-w) - log(t_d[target]-t_a[target]+1) - log(susceptible_population.length()+1) - log(col_pop(data,t_c[target]));
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

List RemoveColonisationTime_List(List data, NumericVector parameters, double &loglik_cur, double w, IntegerVector Vs, IntegerVector &Va, int &nacc_remove)
{
  IntegerVector t_a = data["t_a"];
  IntegerVector t_c = data["t_c"];
  IntegerVector t_d = data["t_d"];
  IntegerVector source = data["source"];
  
  List data_can = clone(data);
  IntegerVector t_c_can = data_can["t_c"];
  IntegerVector source_can = data_can["source"];
  
  IntegerVector susceptible_population = setdiff(Vs,Va);
  IntegerVector nadd0 = ReturnNadd0(data, Va);
  if(nadd0.length()==0) return data;
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
    log_prop_ratio = log(nadd0.length()) + log(1-w) - log(t_d[target]-t_a[target]+1) - log(susceptible_population.length()+1) - log(col_pop(data,t_c[target]));
    //log_prop_ratio = log(Va.length()) + log(1-w) - log(t_d[target]-t_a[target]+1) - log(susceptible_population.length()+1) - log(col_pop(data,t_c[target]));
    
  }
  
  //if(!CheckData(data_can)) stop("error with data in the remove step!");
  bool impute_distances = false;
  if(data.containsElementNamed("genetic_ids") && impute_distances)
  {
    Rcout << "Before update imputed nodes" << std::endl;
    data_can = UpdateImputedNodes(data_can, parameters);
    Rcout << "After update imputed nodes" << std::endl;
    double genetic_contribution = data_can["genetic_contribution"];
    //if(print_output) Rcout << "Gen contribution = " << genetic_contribution << std::endl;
    log_prop_ratio += genetic_contribution;
  }
  
  //Rcout << "Before likelihood" << std::endl;
  double loglik_can = CalculateLogLikelihood(data_can, parameters);
  //Rcout << "After likelihood" << std::endl;
  double U = R::runif(0.0,1.0);
  
  if(log(U) < loglik_can - loglik_cur + log_prop_ratio)
  {
    
    loglik_cur = loglik_can;
    int Va_loc = WhichVec(target, Va)[0];
    Va = RemoveElement(Va_loc, Va);
    //Rcout << " - accepted. Number added = " << Va.length() << std::endl;
    nacc_remove++;
    return data_can;
  }
  else
  {
    return data;
    //Rcout << " - rejected." << std::endl;
  }
}


void ChangeSource(List &data, NumericVector parameters, double &loglik_cur, IntegerVector Vq, IntegerVector Va, int &nacc_source_change)
{
  IntegerVector source = data["source"];
  
  List data_can = clone(data);
  IntegerVector source_can = data_can["source"];
  
  IntegerVector colonised_population = union_(Vq,Va);
  int target = SampleVector(colonised_population);
  if(source[target] == -1) return;
  IntegerVector possible_infectors = ReturnPossibleInfectors(data_can, target);
  if(possible_infectors.length()==0) stop("No infectors!");
  source_can[target] = SampleVector(possible_infectors);
  
  if(source_can[target] == source[target]) return;
  
  double loglik_can = CalculateLogLikelihood(data_can, parameters);
  double U = R::runif(0.0,1.0);
  //Rcout << "Target = " << target << ", previous source = " << source[target] << ", new source = " << source_can[target] << ", loglik = "
  //      << loglik_cur << ", loglik_can = " << loglik_can << std::endl;
  if(log(U) < loglik_can - loglik_cur)
  {
    data = data_can;
    loglik_cur = loglik_can;
    nacc_source_change++;
  }
}

void SwapSourceWithOffspring(List &data, NumericVector parameters, double &loglik_cur, IntegerVector Vq, IntegerVector Va, int &nacc_swap)
{
  IntegerVector t_a = data["t_a"];
  IntegerVector t_d = data["t_d"];
  IntegerVector t_c = data["t_c"];
  IntegerVector source = data["source"];
  
  List data_can = clone(data);
  IntegerVector source_can = data_can["source"];
  IntegerVector t_c_can = data_can["t_c"];
  
  IntegerVector colonised_population = union_(Vq,Va);
  int target = SampleVector(colonised_population);
  if(source[target] == -1) return;
  int target_source = source[target];
  

  IntegerVector target_source_colonisations = WhichVec(target_source,source);
  //Rcout << "Swap target = " << target << " at time " << t_c[target] << " with source = " << target_source << " at time " << t_c[target_source] 
  //      << " source colonisations = " << target_source_colonisations << std::endl;
  for(int i = 0; i<target_source_colonisations.length(); i++)
  {
    int current_source_colonisation = target_source_colonisations[i];
    //Rcout << "current colonisation = " << current_source_colonisation << ", at time t = " << t_c[current_source_colonisation] << std::endl;
    if(t_c[current_source_colonisation] <= t_c[target] && current_source_colonisation != target) return;
  }
  
  int source_coltime = t_c[target_source];
  int target_coltime = t_c[target];
  
  t_c_can[target] = source_coltime;
  t_c_can[target_source] = target_coltime;
  if(t_c_can[target] < t_a[target]) return;
  if(t_c_can[target_source] > t_d[target_source]) return;
  source_can[target] = source[target_source];
  source_can[target_source] = target;
  
  if(source[target] == -1 && t_c[target] != t_a[target])
  {
    stop("import is not infected on day one!");
  }
  
  double loglik_can = CalculateLogLikelihood(data_can, parameters);
  double U = R::runif(0.0,1.0);
  if(log(U) < loglik_can - loglik_cur)
  {
    //Rcout << "Accept the swap!" << std::endl;
    data = data_can;
    loglik_cur = loglik_can;
    nacc_swap++;
  }
  
}

void SwapImportWithSource(List &data, NumericVector parameters, double &loglik_cur, IntegerVector Vq, IntegerVector Va, int &nacc_swap)
{
  IntegerVector t_a = data["t_a"];
  IntegerVector t_c = data["t_c"];
  IntegerVector source = data["source"];
  
  List data_can = clone(data);
  IntegerVector t_c_can = data_can["t_c"];
  IntegerVector source_can = data_can["source"];
  
  IntegerVector importations = WhichVec(-1,source);
  int target = SampleVector(importations);
  IntegerVector target_infections = WhichVec(target,source);
  IntegerVector target_infections_times = as<IntegerVector>(t_c[target_infections]);
  //Rcout << "Target = " << target << std::endl
  //     << "Target infections = " << target_infections << std::endl
  //      << "Target infection times = " << target_infections_times << std::endl;
  if(target_infections.length()==0) return;
  int min_time = target_infections_times[0];
  int min_id = target_infections[0];
  for(int i = 1; i<target_infections.length(); i++)
  {
    if(target_infections_times[i] < min_time)
    {
      min_time = target_infections_times[i];
      min_id = target_infections[i];
    }
  }
  
  // check for duplicates
  IntegerVector min_time_locs = WhichVec(min_time, target_infections_times);
  if(min_time_locs.length() > 1) return;
  
  // update candidate data
  source_can[target] = min_id;
  source_can[min_id] = -1;
  t_c_can[target] = min_time;
  t_c_can[min_id] = t_a[min_id];
  
  if(t_a[target] > min_time) return;
  
  double loglik_can = CalculateLogLikelihood(data_can, parameters);
  double U = R::runif(0.0,1.0);
  //Rcout << "Loglik can = " << loglik_can << ", loglik cur = " << loglik_cur << ", swap import " << target << " with " << min_id << std::endl;
  //Rcout << "t_a[target] = " << t_a[target] << ", t_c[min_id] = " << t_c[min_id] << ", t_a[min_id] = " << t_a[min_id] << std::endl;
  //bool isnan = is_nan();
  if(Rcpp::traits::is_nan<REALSXP>(loglik_can))
  {
    stop("loglik is not a number");
  }
  if(log(U) < loglik_can - loglik_cur)
  {
    //Rcout << "Accept the swap!" << std::endl;
    data = data_can;
    loglik_cur = loglik_can;
    nacc_swap++;
  }
  
}

void ChangeSource_NS(List &data, NumericVector parameters, double &loglik_cur, IntegerVector Vq, IntegerVector Va)
{
  List data_can = clone(data);
  IntegerVector source_can = data_can["source"];
  
  IntegerVector colonised_population = union_(Vq,Va);
  int target = SampleVector(colonised_population);
  
  IntegerVector possible_infectors = ReturnPossibleInfectors(data_can, target);
  if(possible_infectors.length()==0) return;
  source_can[target] = SampleVector(possible_infectors);
  
  double loglik_can = CalculateLogLikelihood_NS(data_can, parameters);
  double U = R::runif(0.0,1.0);
  if(log(U) < loglik_can - loglik_cur)
  {
    data = data_can;
    loglik_cur = loglik_can;
  }
}


void MoveColonisationTime_NS(List &data, NumericVector parameters, double &loglik_cur, double w, IntegerVector Vq, IntegerVector Va, int &nacc_move)
{
  IntegerVector t_a = data["t_a"];
  IntegerVector t_d = data["t_d"];
  IntegerVector t_c = data["t_c"];
  IntegerVector source = data["source"];
  
  List data_can = clone(data);
  IntegerVector t_c_can = data_can["t_c"];
  IntegerVector source_can = data_can["source"];
  
  IntegerVector colonised_population = union_(Vq,Va);
  int target = SampleVector(colonised_population);
  int last_day = t_d[target];
  
  
  
  //Rcout << "Last day = " << last_day << std::endl;
  IntegerMatrix screening_results = data["screening_matrix"];
  IntegerVector row = screening_results(target,_);
  IntegerVector positive_swabs = WhichVec(1,row);
  for(int i = 0; i<positive_swabs.length(); i++)
  {
    int day = positive_swabs[i];
    if(day < last_day)
    {
      last_day = day;
    }
  }
  //Rcout << "Last day = " << last_day << ", positive swabs = " << positive_swabs << std::endl;
  double log_prop_ratio = 0;
  
  double U = R::runif(0.0,1.0);
  if(U < w)
  {
    // Propose the target is an importation
    t_c_can[target] = t_a[target];
    source_can[target] = -1;
    if(source[target] == -1)
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
    t_c_can[target] = SampleVector(seq(t_a[target],last_day));
    IntegerVector possible_infectors = ReturnPossibleInfectors(data_can, target);
    if(possible_infectors.length()==0) return;
    source_can[target] = SampleVector(possible_infectors);
    if(source[target] != -1)
    {
      // Acquisition -> Acquisition
      // do nothing
    }
    else
    {
      // Importation -> Acquisition
      log_prop_ratio = log(w) + log(last_day - t_a[target] + 1) - log(1-w);
    }
  }
  
  if(t_c[target] == t_c_can[target] && source[target] == source_can[target]) return;
  
  
  double loglik_can = CalculateLogLikelihood_NS(data_can, parameters);
  U = R::runif(0.0,1.0);
  //Rcout << "t_a = " << t_a[target] << ", t_d = " << t_d[target] << ", last_day = " << last_day << ", proposed coltime = " << t_c_can[target] << std::endl;
  //Rcout << "  MOVE ID " << target <<  ": Loglik can = " << loglik_can << ", loglik cur = " << loglik_cur << ", log prop ratio = " << log_prop_ratio;
  if(log(U) < loglik_can - loglik_cur + log_prop_ratio)
  {
    data = data_can;
    loglik_cur = loglik_can;
    nacc_move++;
    //Rcout << " - accepted." << std::endl;
  }
  else
  {
    //Rcout << " - rejected." << std::endl;
  }
}

void AddColonisationTime_NS(List &data, NumericVector parameters, double &loglik_cur, double w, IntegerVector Vs, IntegerVector &Va, int &nacc_add)
{
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
    
    log_prop_ratio = log(susceptible_population.length()) - log(w) - log(1+Va.length());
    //log_prop_ratio = log(susceptible_population.length()) - log(w) - log(1+Va.length());
  }
  else
  {
    // Propose an acquisition
    t_c_can[target] = SampleVector(seq(t_a[target],t_d[target]));
    //t_c_can[target] = floor(R::runif(t_a[target],t_d[target]+1));
    IntegerVector possible_infectors = ReturnPossibleInfectors(data_can, target);
    if(possible_infectors.length()==0) return;
    source_can[target] = SampleVector(possible_infectors);
    
    log_prop_ratio = log(susceptible_population.length()) + log(t_d[target] - t_a[target] + 1) - log(1-w) - log(1+Va.length());
  }
  
  
  double loglik_can = CalculateLogLikelihood_NS(data_can, parameters);
  U = R::runif(0.0,1.0);
  //Rcout << "ADDING ID " << target <<  ": Loglik can = " << loglik_can << ", loglik cur = " << loglik_cur << ", log prop ratio = " << log_prop_ratio;
  if(log(U) < loglik_can - loglik_cur + log_prop_ratio)
  {
    data = data_can;
    loglik_cur = loglik_can;
    Va.push_back(target);
    nacc_add++;
    //Rcout << " - accepted. Number added = " << Va.length() << std::endl;
  }
  else
  {
    //Rcout << " - rejected." << std::endl;
  }
}

void RemoveColonisationTime_NS(List &data, NumericVector parameters, double &loglik_cur, double w, IntegerVector Vs, IntegerVector &Va, int &nacc_remove)
{
  IntegerVector t_a = data["t_a"];
  IntegerVector t_c = data["t_c"];
  IntegerVector t_d = data["t_d"];
  IntegerVector source = data["source"];
  
  List data_can = clone(data);
  IntegerVector t_c_can = data_can["t_c"];
  IntegerVector source_can = data_can["source"];
  
  IntegerVector susceptible_population = setdiff(Vs,Va);
  if(Va.length()==0) return;

  int target = SampleVector(Va);
  double log_prop_ratio = 0;
  
  source_can[target] = -2;
  t_c_can[target] = -1;
  
  if(source[target]==-1)
  {
    // Removing an importation
    log_prop_ratio = log(Va.length()) + log(w) - log(susceptible_population.length()+1);
    //log_prop_ratio = log(Va.length()) + log(w) - log(susceptible_population.length()+1);
  }
  else
  {
    // Removing an acquisition
    log_prop_ratio = log(Va.length()) + log(1-w) - log(t_d[target]-t_a[target]+1) - log(susceptible_population.length()+1);
    //log_prop_ratio = log(Va.length()) + log(1-w) - log(t_d[target]-t_a[target]+1) - log(susceptible_population.length()+1) - log(col_pop(data,t_c[target]));
    
  }
  
  
  double loglik_can = CalculateLogLikelihood_NS(data_can, parameters);
  double U = R::runif(0.0,1.0);
  //Rcout << "REMOVE ID " << target <<  ": Loglik can = " << loglik_can << ", loglik cur = " << loglik_cur << ", log prop ratio = " << log_prop_ratio;
  if(log(U) < loglik_can - loglik_cur + log_prop_ratio)
  {
    data = data_can;
    loglik_cur = loglik_can;
    int Va_loc = WhichVec(target, Va)[0];
    Va = RemoveElement(Va_loc, Va);
    nacc_remove++;
    //Rcout << " - accepted. Number added = " << Va.length() << std::endl;
  }
  else
  {
    //Rcout << " - rejected." << std::endl;
  }
}

// [[Rcpp::export]]
void MCMC_Epi(List MCMC_options, IntegerVector t_a, IntegerVector t_c, IntegerVector t_d, IntegerVector source, IntegerMatrix screening_matrix)
{
  // Initialise data
  List data = InitialiseData(t_a, t_c, t_d, source, screening_matrix);
  
  // Load MCMC options
  const int iterations = MCMC_options["iterations"];
  const int population_size = t_a.length();
  const int num_updates = MCMC_options["num_updates"];
  NumericVector prior_parameters = MCMC_options["prior_parameters"];
  NumericVector initial_chain_state = MCMC_options["initial_chain_state"];
  NumericVector parameters = clone(initial_chain_state);
  NumericVector proposal_variance = MCMC_options["proposal_variance"];
  IntegerVector debug_flags = MCMC_options["debug_flags"];
  std::string output_file = MCMC_options["output_file"];
  std::string coltime_file = MCMC_options["coltime_file"];
  std::string source_file = MCMC_options["source_file"];
  std::string loglik_file = MCMC_options["loglik_file"];
  
  // Calculate other variables
  int true_positives = CalculateTruePositives(data);
  int nacc_beta = 0;
  int nacc_add = 0;
  int nacc_move = 0;
  int nacc_remove = 0;
  int nacc_source_change = 0;
  int nacc_swap = 0;
  int nacc_import_swap = 0;
  IntegerVector augmented_moves_proposed(5);
  double loglik = CalculateLogLikelihood(data, parameters);
  IntegerVector Va;
  IntegerVector Vs = CalculateVs(data);
  IntegerVector Vq = CalculateVq(data);
  
  // Write to files
  // Output file
  remove(output_file.c_str());
  std::ofstream myfile; // Define output stream
  myfile.open(output_file.c_str()); // Open file
  assert(myfile.is_open());
  PrintNumVectorToFile(myfile, parameters);
  
  // Colonisation time sum file
  remove(coltime_file.c_str());
  std::ofstream myfile2;
  myfile2.open(coltime_file.c_str());
  assert(myfile2.is_open());
  PrintColonisationTimeSumToFile(myfile2, data);
  
  // Source file
  remove(source_file.c_str());
  std::ofstream myfile3;
  myfile3.open(source_file.c_str());
  assert(myfile3.is_open());
  PrintIntVectorToFile(myfile3, data["source"]);
  
  // loglik file
  remove(loglik_file.c_str());
  std::ofstream myfile4;
  myfile4.open(loglik_file.c_str());
  assert(myfile4.is_open());
  myfile4 << loglik << std::endl;
  
  
  // Begin MCMC
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
      int import_sum = CalculateImportSum(data);
      parameters[1] = R::rbeta(prior_parameters[2] + import_sum, prior_parameters[3] + population_size - import_sum);
    }
    
    loglik = CalculateLogLikelihood(data, parameters);
    
    // Update beta by metropolis hastings
    if(debug_flags[2]==0)
    {
      UpdateTransmissionRate(data, parameters, loglik, proposal_variance[0], prior_parameters[4], nacc_beta);
    }
    
    // Augmented data updates
    if(debug_flags[4]==0)
    {
      for(int update_counter = 0; update_counter < num_updates; update_counter++)
      {
        double w = 0.3;
        double move = R::runif(0.0,1.0);
        augmented_moves_proposed[floor(move)]++;
        if(floor(move) < 1)
        {
          // Move a colonisation time
          MoveColonisationTime(data, parameters, loglik, w, Vq, Va, nacc_move);
        }
        else if(floor(move) < 2)
        {
          // Add a colonsation time
          AddColonisationTime(data, parameters, loglik, w, Vs, Va, nacc_add);
        }
        else if(floor(move) < 3)
        {
          // Remove a colonisation time
          RemoveColonisationTime(data, parameters, loglik, w, Vs, Va, nacc_remove);
        }
        else if(floor(move) < 4)
        {
          // change a source without changing a colonisation time
          ChangeSource(data, parameters, loglik, Vq, Va, nacc_source_change);
        }
        else if(floor(move) < 5)
        {
          // swap a source with the offspring
          SwapSourceWithOffspring(data, parameters, loglik, Vq, Va, nacc_swap);
        }
        else if(floor(move) < 6)
        {
          // swap an import with the first possible source
          SwapImportWithSource(data, parameters, loglik, Vq, Va, nacc_import_swap);
        }
        
        // Write colonisation time sum to file
        PrintColonisationTimeSumToFile(myfile2, data);
        PrintIntVectorToFile(myfile3, data["source"]);
        
        myfile4 << loglik << std::endl;
      }
    }
      
    
    
    // Write parameters to file
    PrintNumVectorToFile(myfile, parameters);
    
    if(i%100==0)
    {
      Rcout << i << std::endl;
    }
    
  }
  
  // Close files
  myfile.close();
  myfile2.close();
  myfile3.close();
  myfile4.close();
  
  double prob = (double)nacc_beta/(double)iterations;
  double move_prob = (double)nacc_move/(double)(augmented_moves_proposed[0]);
  double add_prob = (double)nacc_add/(double)(augmented_moves_proposed[1]);
  double remove_prob = (double)nacc_remove/(double)(augmented_moves_proposed[2]);
  double source_change_prob = (double)nacc_source_change/(double)(augmented_moves_proposed[3]);
  double swap_prob = (double)nacc_swap/(double)(augmented_moves_proposed[4]);
  double import_swap_prob = (double)nacc_import_swap/(double)(augmented_moves_proposed[5]);
  Rcout << "Beta acceptance = " << prob << std::endl;
  Rcout << "Move prob = " << move_prob << ", add prob = " << add_prob << ", remove prob = " << remove_prob 
        << ", source change prob = " << source_change_prob << ", swap prob = " << swap_prob 
        << ", import swap prob = " << import_swap_prob << std::endl;
  
  
  
}


void PrintMatrix(IntegerMatrix matrix)
{
  for(int i = 0; i<matrix.nrow(); i++)
  {
    for(int j = 0; j<matrix.ncol(); j++)
    {
      Rcout << matrix(i,j) << " ";
    }
    Rcout << std::endl;
  }
}



// [[Rcpp::export]]
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


bool DoesTargetHaveSwabAtTime(List data, int target, int current_variant, int time)
{
  IntegerVector genetic_ids = data["genetic_ids"];
  IntegerVector variant_numbers = data["variant_numbers"];
  IntegerVector sample_times = data["sample_times"];
  IntegerVector t_c = data["t_c"];
  
  IntegerVector target_swab_loc = WhichVec(target, genetic_ids);
  IntegerVector variant_loc = WhichVec(current_variant, variant_numbers);
  IntegerVector target_variant_swab_locs = intersect(target_swab_loc, variant_loc);
  IntegerVector target_swab_times = as<IntegerVector>(sample_times[target_variant_swab_locs]);
  IntegerVector target_swab_loc_same_time = WhichVec(time, target_swab_times);
  
  if(target_swab_loc_same_time.length()==0)
  {
    // There is no swab at the time
    return false;
  }
  else
  {
    // There is a swab at the time
    return true;
  }
}



// This function will take the current configuration of the transmission tree, and sequentially add nodes that need to be imputed
List ImputeNodes(List data, NumericVector parameters)
{
  bool print_debug = true;
  //List data_can = clone(data);
  IntegerVector genetic_ids = data["genetic_ids"];
  IntegerVector sample_times = data["sample_times"];
  IntegerMatrix genetic_matrix = data["genetic_matrix"];
  IntegerVector variant_numbers = data["variant_numbers"];
  
  int original_length = genetic_ids.length();
  
  IntegerVector source = data["source"];
  IntegerVector t_c = data["t_c"];
  
  double mutation_rate = parameters[3];
  
  IntegerVector ever_infected;
  if(print_debug)
  {
    Rcout << "Try to impute nodes" << std::endl;
    Rcout << "Genetic ids  = " << genetic_ids << std::endl;
    Rcout << "Sample times = " << sample_times << std::endl;
  }
  
  
  for(int i = 0; i<t_c.length(); i++)
  {
    if(t_c[i] != -1 && source[i] != -1)
    {
      // Check if they have a swab at the time of infection or if their infector has a swab or another infection
      int patient_id = i;
      int current_variant = 1;
      int col_time = t_c[i];
      if(print_debug) Rcout << "Check patient " << patient_id << ", with t_c = " << col_time << std::endl;
      IntegerVector patient_swab_loc = WhichVec(patient_id, genetic_ids);
      IntegerVector variant_loc = WhichVec(current_variant, variant_numbers);
      IntegerVector patient_variant_swab_locs = intersect(patient_swab_loc, variant_loc);
      IntegerVector patient_swab_times = as<IntegerVector>(sample_times[patient_variant_swab_locs]);
      IntegerVector patient_swab_loc_same_time = WhichVec(col_time, patient_swab_times);
      if(patient_swab_loc_same_time.length()==0)
      {
        if(print_debug) Rcout << "No swab found at the time of infection - checking infector swabs" << std::endl;
        // Now check infector swabs and infector other infections
        int infector_id = source[patient_id];
        IntegerVector infector_swab_loc = WhichVec(infector_id, genetic_ids);
        IntegerVector infector_variant_swab_locs = intersect(infector_swab_loc, variant_loc);
        IntegerVector infector_swab_times = as<IntegerVector>(sample_times[infector_variant_swab_locs]);
        IntegerVector infector_swab_loc_same_time = WhichVec(col_time, infector_swab_times);
        if(infector_swab_loc_same_time.length()==0)
        {
          if(print_debug) Rcout << "No swab found in the infector, checking their infections for swabs at the time of colonisation" << std::endl;
          // The infector also does not have a swab at the time of infection, check if there are any other infections at this time
          IntegerVector infectees = WhichVec(infector_id, source);
          
          bool swab_found = false;
          for(int j = 0; j<infectees.length(); j++)
          {
            int current_infectee = infectees[j];
            if(current_infectee != i && t_c[current_infectee] == t_c[i])
            {
              // The infectee is not the original target, and also has a colonisation at the same time
              IntegerVector infectee_swab_loc = WhichVec(current_infectee, genetic_ids);
              IntegerVector infectee_variant_swab_locs = intersect(infectee_swab_loc, variant_loc);
              IntegerVector infectee_swab_times = as<IntegerVector>(sample_times[infectee_variant_swab_locs]);
              IntegerVector infectee_swab_loc_same_time = WhichVec(col_time, infectee_swab_times);
              if(infectee_swab_loc_same_time.length()!=0)
              {
                swab_found = true;
              }
            }
          }
          
          if(!swab_found)
          {
            if(print_debug) Rcout << "No swab has been found, therefore we must impute distances" << std::endl;
            // there has been no swab found at the same time, therefore try and add a node
            IntegerVector updated_genetic_ids = clone(genetic_ids);
            IntegerVector updated_sample_times = clone(sample_times);
            IntegerVector updated_variant_numbers = clone(variant_numbers);
            updated_genetic_ids.push_back(i);
            updated_sample_times.push_back(col_time);
            updated_variant_numbers.push_back(1);
            List temp_data = List::create(Named("t_c") = t_c,
                                          Named("source") = source,
                                          Named("genetic_ids") = updated_genetic_ids,
                                          Named("sample_times") = updated_sample_times,
                                          Named("variant_numbers") = updated_variant_numbers);
            IntegerVector updated_gen_source = ReturnGenSourceVector_WHD(temp_data);
            if(print_debug) Rcout << "Updated genetic source vector = " << updated_gen_source << std::endl;
            IntegerVector children_nodes = WhichVec(updated_gen_source.length()-1, updated_gen_source);
            if(print_debug) Rcout << "Children nodes = " << children_nodes << std::endl;
            if(children_nodes.length() > 1)
            {
              // The imputed node is valid (has 2 or more links, therefore impute matrix)
              int parent_node = updated_gen_source[updated_gen_source.length()-1];
              
              if(parent_node == -1)
              {
                // We are imputing an imported case
                // Find the child node with the minimum time
                int min_node = children_nodes[0];
                int min_time = sample_times[min_node];
                for(int j = 1; j<children_nodes.length(); j++)
                {
                  int current_child = children_nodes[j];
                  int current_time = sample_times[current_child];
                  if(current_time < min_time)
                  {
                    min_node = current_child;
                    min_time = current_time;
                  }
                }
                double min_time_double = min_time;
                int draw = CalculateMutationProbability(mutation_rate, min_time_double); 
                
                // initialise a new genetic matrix of zeros
                IntegerMatrix updated_genetic_matrix(updated_gen_source.length());
                
                // copy entries from the old matrix
                for(int ii = 0; ii<genetic_ids.length(); ii++)
                {
                  for(int jj = 0; jj<genetic_ids.length(); jj++)
                  {
                    updated_genetic_matrix(ii,jj) = genetic_matrix(ii,jj);
                  }
                }
                
                int imputed_node = updated_genetic_ids.length() - 1;
                if(print_debug) Rcout << "Imputed node = " << imputed_node << std::endl;
                // now fill in final column
                updated_genetic_matrix(min_node,imputed_node) = draw;
                updated_genetic_matrix(imputed_node,min_node) = draw;
                
                
                // now fill in the rest
                for(int j = 0; j<updated_gen_source.length()-1; j++)
                {
                  if(j != min_node)
                  {
                    IntegerVector path = ReturnPathToRoot_Rcpp(updated_gen_source, j);
                    if(print_debug)
                    {
                      Rcout << "Node = " << j << ", path = " << path << std::endl;
                    }
                    bool path_contains_min_node = false;
                    for(int k = 0; k<path.length(); k++)
                    {
                      if(path[k] == min_node)
                      {
                        path_contains_min_node = true;
                        break;
                      }
                    }
                    if(path_contains_min_node)
                    {
                      updated_genetic_matrix(j,imputed_node) = updated_genetic_matrix(min_node,j) + updated_genetic_matrix(min_node,imputed_node);
                      updated_genetic_matrix(imputed_node,j) = updated_genetic_matrix(min_node,j) + updated_genetic_matrix(min_node,imputed_node);
                      if(updated_genetic_matrix(j,imputed_node) < 0) stop("-ve distance");
                    }
                    else
                    {
                      updated_genetic_matrix(j,imputed_node) = updated_genetic_matrix(min_node,j) - updated_genetic_matrix(min_node,imputed_node);
                      updated_genetic_matrix(imputed_node,j) = updated_genetic_matrix(min_node,j) - updated_genetic_matrix(min_node,imputed_node);
                      if(updated_genetic_matrix(j,imputed_node) < 0) stop("-ve distance");
                    }
                  }
                }
                if(print_debug)
                {
                  Rcout << "Printing updated genetic matrix" << std::endl;
                  for(int ii = 0; ii<updated_genetic_matrix.nrow(); ii++)
                  {
                    for(int jj = 0; jj<updated_genetic_matrix.ncol(); jj++)
                    {
                      Rcout << updated_genetic_matrix(ii,jj) << " ";
                    }
                    Rcout << std::endl;
                  }
                }
                
                
                genetic_ids = updated_genetic_ids;
                sample_times = updated_sample_times;
                genetic_matrix = updated_genetic_matrix;
                variant_numbers = updated_variant_numbers;
                
              }
              else
              {
                // We are importing an internal node
                int min_node = children_nodes[0];
                int min_distance = genetic_matrix(parent_node, children_nodes[0]);
                for(int j = 1; j<children_nodes.length(); j++)
                {
                  int current_child = children_nodes[j];
                  int current_distance = genetic_matrix(parent_node, current_child);
                  if(current_distance < min_distance)
                  {
                    min_node = current_child;
                    min_distance = current_distance;
                  }
                }
                if(print_debug)
                {
                  Rcout << "Parent node = " << parent_node << ", min node = " << min_node << ", distance = " << min_distance << std::endl;
                }
                
                int max_time_diff = sample_times[min_node] - sample_times[parent_node];
                int current_time_diff = col_time - sample_times[parent_node];
                double time_ratio = (double)(current_time_diff)/(double)(max_time_diff);
                if(max_time_diff == 0) stop("error in logic, there is a node at the same time");
                
                int draw = R::rbinom(min_distance, time_ratio);
                if(print_debug) Rcout << "max time diff = " << max_time_diff << ", current time diff = " 
                                      << current_time_diff << ", ratio = " << time_ratio << std::endl;
                // now add a new column in the genetic matrix
                
                // initialise a new genetic matrix of zeros
                IntegerMatrix updated_genetic_matrix(updated_gen_source.length());
                
                // copy entries from the old matrix
                for(int ii = 0; ii<genetic_ids.length(); ii++)
                {
                  for(int jj = 0; jj<genetic_ids.length(); jj++)
                  {
                    updated_genetic_matrix(ii,jj) = genetic_matrix(ii,jj);
                  }
                }
                
                int imputed_node = updated_genetic_ids.length() - 1;
                if(print_debug) Rcout << "Imputed node = " << imputed_node << std::endl;
                // now fill in final column
                updated_genetic_matrix(parent_node,imputed_node) = draw;
                updated_genetic_matrix(imputed_node,parent_node) = draw;
                
                // now fill in the rest
                for(int j = 0; j<updated_gen_source.length()-1; j++)
                {
                  if(j != parent_node)
                  {
                    IntegerVector path = ReturnPathToRoot_Rcpp(updated_gen_source, j);
                    if(print_debug) Rcout << "Node = " << j << ", path = " << path << std::endl;
                    bool path_contains_imputed_node = false;
                    for(int k = 0; k<path.length(); k++)
                    {
                      if(path[k] == imputed_node)
                      {
                        path_contains_imputed_node = true;
                        break;
                      }
                    }
                    if(path_contains_imputed_node)
                    {
                      updated_genetic_matrix(j,imputed_node) = updated_genetic_matrix(parent_node,j) - updated_genetic_matrix(parent_node,imputed_node);
                      updated_genetic_matrix(imputed_node,j) = updated_genetic_matrix(parent_node,j) - updated_genetic_matrix(parent_node,imputed_node);
                      if(updated_genetic_matrix(j,imputed_node) < 0)
                      {
                        Rcout << "Parent node = " << parent_node << ", imputed node = " << imputed_node << std::endl;
                        Rcout << "Dist 1 = " << updated_genetic_matrix(parent_node,j) << ", dist 2 = " << updated_genetic_matrix(parent_node,imputed_node) << std::endl;
                        stop("-ve distance");
                      }
                    }
                    else
                    {
                      updated_genetic_matrix(j,imputed_node) = updated_genetic_matrix(parent_node,j) + updated_genetic_matrix(parent_node,imputed_node);
                      updated_genetic_matrix(imputed_node,j) = updated_genetic_matrix(parent_node,j) + updated_genetic_matrix(parent_node,imputed_node);
                    }
                  }
                }
                if(print_debug)
                {
                  Rcout << "Printing updated genetic matrix" << std::endl;
                  for(int ii = 0; ii<updated_genetic_matrix.nrow(); ii++)
                  {
                    for(int jj = 0; jj<updated_genetic_matrix.ncol(); jj++)
                    {
                      Rcout << updated_genetic_matrix(ii,jj) << " ";
                    }
                    Rcout << std::endl;
                  }
                }
                
                
                genetic_ids = updated_genetic_ids;
                sample_times = updated_sample_times;
                genetic_matrix = updated_genetic_matrix;
                variant_numbers = updated_variant_numbers;
              }
            }
          }
        }
      }
      
    }
  }
  
  int new_length = genetic_ids.length();
  IntegerVector imputed_nodes(new_length);
  for(int i = original_length; i<new_length; i++)
  {
    imputed_nodes[i] = 1;
  }
  
  if(print_debug) Rcout << "End of fn - genetic ids = " << genetic_ids << std::endl;
  List return_list = List::create(Named("t_a") = data["t_a"], 
                                  Named("t_c") = t_c,
                                  Named("t_d") = data["t_d"],
                                  Named("source") = source,
                                  Named("screening_matrix") = data["screening_matrix"],
                                  Named("genetic_ids") = genetic_ids,
                                  Named("sample_times") = sample_times,
                                  Named("variant_numbers") = variant_numbers,
                                  Named("genetic_matrix") = genetic_matrix,
                                  Named("sequence_length") = data["sequence_length"],
                                  Named("imputed_nodes") = imputed_nodes);
  return return_list;
}

// Given data already containing some imputed nodes, we wish to update the genetic side given different epi data
// and also return the proposal ratio probability
List UpdateImputedNodes_old(List data, NumericVector parameters)
{
  bool print_debug = false;
  // First we wish to use ONLY the observed data, and work out nodes that need to be imputed
  IntegerVector t_c = data["t_c"];
  IntegerVector source = data["source"];
  IntegerVector sample_times = data["sample_times"];
  IntegerVector imputed_nodes = data["imputed_nodes"];
  IntegerVector genetic_ids = data["genetic_ids"];
  IntegerVector variant_numbers = data["variant_numbers"];
  IntegerMatrix genetic_matrix = data["genetic_matrix"];
  
  IntegerVector observed_loc = WhichVec(0,imputed_nodes);
  IntegerVector observed_genetic_ids = as<IntegerVector>(genetic_ids[observed_loc]);
  IntegerVector observed_sample_times = as<IntegerVector>(sample_times[observed_loc]);
  IntegerVector observed_variant_numbers = as<IntegerVector>(variant_numbers[observed_loc]);
  IntegerMatrix observed_genetic_matrix(observed_loc.length());
  
  IntegerVector genetic_ids_can = clone(observed_genetic_ids);
  IntegerVector sample_times_can = clone(observed_sample_times);
  IntegerVector variant_numbers_can = clone(observed_variant_numbers);
  IntegerMatrix genetic_matrix_can = clone(observed_genetic_matrix);
  
  for(int i = 0; i<observed_loc.length(); i++)
  {
    for(int j = 0; j<observed_loc.length(); j++)
    {
      int ii = observed_loc[i];
      int jj = observed_loc[j];
      observed_genetic_matrix(i,j) = genetic_matrix(ii,jj);
    }
  }
  
  // Now return a list of updated genetic ids
  
  for(int i = 0; i<t_c.length(); i++)
  {
    if(t_c[i] != -1 && source[i] != -1)
    {
      // Check if they have a swab at the time of infection or if their infector has a swab or another infection
      int patient_id = i;
      int current_variant = 1;
      int col_time = t_c[i];
      if(print_debug) Rcout << "Check patient " << patient_id << ", with t_c = " << col_time << std::endl;
      IntegerVector patient_swab_loc = WhichVec(patient_id, observed_genetic_ids);
      IntegerVector variant_loc = WhichVec(current_variant, observed_variant_numbers);
      IntegerVector patient_variant_swab_locs = intersect(patient_swab_loc, variant_loc);
      IntegerVector patient_swab_times = as<IntegerVector>(observed_sample_times[patient_variant_swab_locs]);
      IntegerVector patient_swab_loc_same_time = WhichVec(col_time, patient_swab_times);
      if(patient_swab_loc_same_time.length()==0)
      {
        if(print_debug) Rcout << "No swab found at the time of infection - checking infector swabs" << std::endl;
        // Now check infector swabs and infector other infections
        int infector_id = source[patient_id];
        IntegerVector infector_swab_loc = WhichVec(infector_id, observed_genetic_ids);
        IntegerVector infector_variant_swab_locs = intersect(infector_swab_loc, variant_loc);
        IntegerVector infector_swab_times = as<IntegerVector>(observed_sample_times[infector_variant_swab_locs]);
        IntegerVector infector_swab_loc_same_time = WhichVec(col_time, infector_swab_times);
        if(infector_swab_loc_same_time.length()==0)
        {
          if(print_debug) Rcout << "No swab found in the infector, checking their infections for swabs at the time of colonisation" << std::endl;
          // The infector also does not have a swab at the time of infection, check if there are any other infections at this time
          IntegerVector infectees = WhichVec(infector_id, source);
          
          bool swab_found = false;
          for(int j = 0; j<infectees.length(); j++)
          {
            int current_infectee = infectees[j];
            if(current_infectee != i && t_c[current_infectee] == t_c[i])
            {
              // The infectee is not the original target, and also has a colonisation at the same time
              IntegerVector infectee_swab_loc = WhichVec(current_infectee, observed_genetic_ids);
              IntegerVector infectee_variant_swab_locs = intersect(infectee_swab_loc, variant_loc);
              IntegerVector infectee_swab_times = as<IntegerVector>(observed_sample_times[infectee_variant_swab_locs]);
              IntegerVector infectee_swab_loc_same_time = WhichVec(col_time, infectee_swab_times);
              if(infectee_swab_loc_same_time.length()!=0)
              {
                swab_found = true;
              }
            }
          }
          
          if(!swab_found)
          {
            if(print_debug) Rcout << "No swab has been found, therefore we must impute distances" << std::endl;
            // there has been no swab found at the same time, therefore try and add a node
            IntegerVector updated_genetic_ids = clone(genetic_ids_can);
            IntegerVector updated_sample_times = clone(sample_times_can);
            IntegerVector updated_variant_numbers = clone(variant_numbers_can);
            updated_genetic_ids.push_back(i);
            updated_sample_times.push_back(col_time);
            updated_variant_numbers.push_back(1);
            List temp_data = List::create(Named("t_c") = t_c,
                                          Named("source") = source,
                                          Named("genetic_ids") = updated_genetic_ids,
                                          Named("sample_times") = updated_sample_times,
                                          Named("variant_numbers") = updated_variant_numbers);
            IntegerVector updated_gen_source = ReturnGenSourceVector_WHD(temp_data);
            if(print_debug) Rcout << "Updated genetic source vector = " << updated_gen_source << std::endl;
            IntegerVector children_nodes = WhichVec(updated_gen_source.length()-1, updated_gen_source);
            if(print_debug) Rcout << "Children nodes = " << children_nodes << std::endl;
            if(children_nodes.length() > 1)
            {
              // the imputed node is valid, so add it to the vector of genetic ids and sample_times
              genetic_ids_can = updated_genetic_ids;
              sample_times_can = updated_sample_times;
              variant_numbers_can = updated_variant_numbers;
            }
          }
        }
      }
      
    }
  }
  
  List temp_data_can = List::create(Named("t_c") = t_c,
                                    Named("source") = source,
                                    Named("genetic_ids") = genetic_ids_can,
                                    Named("sample_times") = sample_times_can,
                                    Named("variant_numbers") = variant_numbers_can);
  // Now compare the links between the original genetic ids and the new ones
  IntegerVector gen_source_can = ReturnGenSourceVector_WHD(temp_data_can);
  IntegerVector gen_source_cur = ReturnGenSourceVector_WHD(data);
  Rcout << "Observed genetic ids = " << observed_genetic_ids << std::endl;
  Rcout << "Genetic ids cur      = " << genetic_ids << std::endl;
  Rcout << "genetic ids can      = " << genetic_ids_can << std::endl;
  Rcout << "Gen source cur       = " << gen_source_cur << std::endl;
  Rcout << "Gen source can       = " << gen_source_can << std::endl;
  
  
  for(int i = observed_loc.length(); i<gen_source_can.length(); i++)
  {
    int parent_node_can = gen_source_can[i];
    int imputed_node_can = i;
    int draw;
    // Check if this combination appears in the current data
    bool already_seen_link = false;
    for(int j = observed_loc.length(); j<gen_source_cur.length(); j++)
    {
      int parent_node_cur = gen_source_cur[j];
      int imputed_node_cur = j;
      if(parent_node_can == parent_node_cur && imputed_node_can == imputed_node_cur)
      {
        already_seen_link = true;
        break;
      }
    }
    
    if(already_seen_link)
    {
      // The node has previous been seen in the data, therefore no need to randomly generate a distance, copy from previous matrix
      if(parent_node_can != -1)
      {
        // the node is internal (not imputed)
        draw = genetic_matrix(parent_node_can, imputed_node_can);
      }
      else
      {
        // the node is imputed
        stop("fix adding imputed nodes");
      }
    }
    else
    {
      // The node has not been seen before, therefore continue as usual and randomly generate distances which are consistant with the minimum
      IntegerVector children_nodes = WhichVec(imputed_node_can, gen_source_can);
      int min_node = children_nodes[0];
      int min_distance = genetic_matrix(parent_node_can, children_nodes[0]);
      for(int j = 1; j<children_nodes.length(); j++)
      {
        int current_child = children_nodes[j];
        int current_distance = genetic_matrix(parent_node_can, current_child);
        if(current_distance < min_distance)
        {
          min_node = current_child;
          min_distance = current_distance;
        }
      }
      if(print_debug)
      {
        Rcout << "Parent node = " << parent_node_can << ", min node = " << min_node << ", distance = " << min_distance << std::endl;
      }
      
      int max_time_diff = sample_times_can[min_node] - sample_times_can[parent_node_can];
      int current_time_diff = sample_times_can[imputed_node_can] - sample_times_can[parent_node_can];
      double time_ratio = (double)(current_time_diff)/(double)(max_time_diff);
      if(max_time_diff == 0) stop("error in logic, there is a node at the same time");
      
      draw = R::rbinom(min_distance, time_ratio);
    }
    
    // now fill in the rest of the matrix deterministically
    
    if(draw==1000) Rcout << "do nothing" << std::endl; // DO NOTHING ---- DELETE THIS
  }
  
  
  
  if(gen_source_can.length() != gen_source_cur.length()) stop("stop function");
  return data;
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

/*

List InitialiseMCMCImputedNodes(List data)
{
  IntegerVector t_c = data["t_c"];
  IntegerVector source = data["source"];
  IntegerVector sample_times = data["sample_times"];
  IntegerVector genetic_ids = data["genetic_ids"];
  IntegerVector variant_numbers = data["variant_numbers"];
  IntegerMatrix original_genetic_matrix = data["genetic_matrix"];
  IntegerMatrix genetic_matrix = clone(original_genetic_matrix);
  int num_patients = t_c.length();
  int current_variant = 1;
  
  
  IntegerVector function_imputed_nodes;
  IntegerVector gen_source_imputed_nodes;
  for(int i = 0; i<num_patients; i++)
  {
    if(source[i] != -1)
    {
      bool does_node_need_impute = DoesNodeNeedToBeImputed(data, i);
      if(does_node_need_impute)
      {
        function_imputed_nodes.push_back(i);
      }
    }

  }

  // Loop through everyone ever infected and check if they need a sequence at the time of colonisation
  for(int i = 0; i<num_patients; i++)
  {
    int target = i;
    int col_time = t_c[target];
    //Rcout << "Target = " << target << ", t_c = " << col_time << std::endl;
    if(t_c[i] != -1 && source[i] != -1)
    {
      // The current target has been colonised
      bool does_target_have_swab_at_col_time = CheckIfThereIsASwabAtTime(data, target, col_time, current_variant);
      if(!does_target_have_swab_at_col_time)
      {
        Rcout << "No time found for target i = " << target << std::endl;
        // There was no swab found at the time of colonisation, therefore check if a node needs to be imputed
        // First we add the genetic ID of the SOURCE (!!) and check the corresponding genetic network, if there are
        // two children nodes, then the node needs to be kept there
        int current_source = source[target];
        Rcout << "current source = " << current_source << std::endl;
        
        IntegerVector updated_genetic_ids = clone(genetic_ids);
        IntegerVector updated_sample_times = clone(sample_times);
        IntegerVector updated_variant_numbers = clone(variant_numbers);
        updated_genetic_ids.push_back(current_source);
        updated_sample_times.push_back(col_time);
        updated_variant_numbers.push_back(current_variant);
        
        List temp_data = List::create(Named("t_c") = t_c,
                                      Named("source") = source,
                                      Named("genetic_ids") = updated_genetic_ids,
                                      Named("sample_times") = updated_sample_times,
                                      Named("variant_numbers") = updated_variant_numbers);
        IntegerVector updated_gen_source = ReturnGenSourceVector_WHD(temp_data);
        int genetic_ids_length = updated_gen_source.length();
        IntegerVector test_node_locs = WhichVec(genetic_ids_length-1, updated_gen_source);
        Rcout << "Gen source = " << updated_gen_source << std::endl;
        Rcout << "Test node locs = " << test_node_locs << std::endl;
        
        if(test_node_locs.length()>1)
        {
          gen_source_imputed_nodes.push_back(i);
          genetic_ids = updated_genetic_ids;
          sample_times = updated_sample_times;
          variant_numbers = updated_variant_numbers;
          
          // Now update the genetic matrix and copy the source sequence
          IntegerMatrix updated_genetic_matrix(genetic_ids_length,genetic_ids_length);
          for(int ii = 0; ii<(genetic_ids_length-1); ii++)
          {
            for(int jj = 0; jj<(genetic_ids_length-1); jj++)
            {
              updated_genetic_matrix(ii,jj) = genetic_matrix(ii,jj);
            }
          }
          
          int imputed_node = genetic_ids_length-1;
          int current_gen_source = updated_gen_source[imputed_node];
          Rcout << "Imputed node = " << imputed_node << ", current gen source " << current_gen_source << std::endl;
          
          if(current_gen_source == -1)
          {
            // find the smallest time and set to zero for initialisation and fill in all other links
            IntegerVector children_nodes = WhichVec(imputed_node, updated_gen_source);
            Rcout << "Children nodes = " << children_nodes << std::endl;
            int min_node = children_nodes[0];
            int min_time = sample_times[min_node] - col_time;
            if(min_time <= 0) stop("min time <= 0");
            for(int j = 0; j<children_nodes.length(); j++)
            {
              int current_node = children_nodes[j];
              int current_node_time = sample_times[current_node] - col_time;
              if(current_node_time < min_time)
              {
                min_time = current_node_time;
                min_node = current_node;
              }
            }
            
            updated_genetic_matrix(min_node, imputed_node) = 0;
            updated_genetic_matrix(imputed_node, min_node) = 0;
            
            for(int j = 0; j<(imputed_node-1);j++)
            {
              updated_genetic_matrix(j, imputed_node) = updated_genetic_matrix(j, min_node);
              updated_genetic_matrix(imputed_node, j) = updated_genetic_matrix(min_node, j);
            }
            
          }
          else
          {
            // otherwise just copy
            for(int j = 0; j<(genetic_ids_length-1); j++)
            {
              updated_genetic_matrix(j,genetic_ids_length-1) = genetic_matrix(j,current_gen_source);
              updated_genetic_matrix(genetic_ids_length-1,j) = genetic_matrix(current_gen_source,j);
              //Rcout << "genetic_matrix(j,current_gen_source) = " << genetic_matrix(j,current_gen_source) << std::endl;
            }
          }
          Rcout << "Current gen source = " << current_gen_source << std::endl;
          Rcout << "genetic ids length = " << genetic_ids_length << std::endl;
          Rcout << "sample times length = " << sample_times.length() << std::endl;

          

          genetic_matrix = updated_genetic_matrix;
          Rcout << "Genetic ids = " << genetic_ids << std::endl;
          Rcout << "Sample times = " << sample_times << std::endl;

          PrintMatrix(genetic_matrix);
        }
        
        
      }
    }


  }
  List other_function_imp_nodes = ReturnNodesToImpute(data);
  IntegerVector other_fn_nodes = other_function_imp_nodes["nodes"];
  IntegerVector other_fn_node_times = other_function_imp_nodes["times"];
  Rcout << "gen source imputed nodes" << gen_source_imputed_nodes << std::endl;
  Rcout << "fucntion imputed nodes" << function_imputed_nodes << std::endl;
  Rcout << "other fn imputed nodes " << other_fn_nodes << std::endl;
  Rcout << "other fn node times " << other_fn_node_times << std::endl;
  if(gen_source_imputed_nodes.length() != function_imputed_nodes.length())
  {
    // could add a check at the time
    //stop("issue with determining which nodes need imputation");
  }
  List temp_data = List::create(Named("t_c") = t_c,
                                Named("source") = source,
                                Named("genetic_ids") = genetic_ids,
                                Named("sample_times") = sample_times,
                                Named("variant_numbers") = variant_numbers);
  IntegerVector updated_gen_source = ReturnGenSourceVector_WHD(temp_data);
  Rcout << "before gen source = " << updated_gen_source << std::endl;
  
  IntegerVector imputed_nodes(genetic_ids.length());
  int original_length = original_genetic_matrix.nrow();
  for(int i = original_length; i<genetic_ids.length(); i++)
  {
    imputed_nodes[i] = 1;
  }
  
  IntegerVector imputed_idx = WhichVec(1,imputed_nodes);
  IntegerVector nodes_to_remove;
  for(int i = 0; i<imputed_idx.length(); i++)
  {
    int current_idx = imputed_idx[i];
    IntegerVector children_nodes = WhichVec(current_idx, updated_gen_source);
    if(children_nodes.length() < 2)
    {
      nodes_to_remove.push_back(current_idx);
    }
  }
  
  Rcout << "Nodes to remove = " << nodes_to_remove << std::endl;
  
  for(int i = 0; i<nodes_to_remove.length(); i++)
  {
    int idx_to_remove = nodes_to_remove[i];
    genetic_ids = RemoveElement(idx_to_remove-i, genetic_ids);
    sample_times = RemoveElement(idx_to_remove-i, sample_times);
    imputed_nodes = RemoveElement(idx_to_remove-i, imputed_nodes);
    variant_numbers = RemoveElement(idx_to_remove-i, variant_numbers);
    genetic_matrix = crossout(genetic_matrix,idx_to_remove-i);
  }
  
  
  
  
  
  IntegerVector t_a = data["t_a"];
  IntegerVector t_d = data["t_d"];
  IntegerMatrix screening_matrix = data["screening_matrix"];
  int sequence_length = data["sequence_length"];
  double importation_p = data["importation_p"];
  IntegerVector master_distances(imputed_nodes.length());
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
                                 Named("importation_p") = importation_p,
                                 Named("master_distances") = master_distances);
  IntegerVector new_gen_source = ReturnGenSourceVector_WHD(final_data);
  
  IntegerVector import_idx = WhichVec(-1, new_gen_source);
  
  for(int i = 0; i<import_idx.length(); i++)
  {
    int current_import_idx = import_idx[i];
    int master_distance = R::rbinom(sequence_length, importation_p);
    master_distances[current_import_idx] = master_distance;
  }
  

  bool print_debug = true;
  if(print_debug)
  {
    Rcout << "Final gen source = " << new_gen_source << std::endl;
  }
  
  
  

  return final_data;
}

*/


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
    IntegerVector gen_source = ReturnGenSourceVector_WHD(temp_data);
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
                                 Named("imputed_nodes") = imputed_nodes);
  IntegerVector new_gen_source = ReturnGenSourceVector_WHD(final_data);
  
  

  

  
  bool print_debug = true;
  if(print_debug)
  {
    Rcout << "Final gen source = " << new_gen_source << std::endl;
  }
 //stop("stop early");
  
  return final_data;
}

//bool DoesTargetNeedANodeImputed(List data, int target, int time, int current_variant)
//{
//}

bool CheckIfThereIsASwabAtTime(List data, int target, int time, int current_variant)
{
  // The idea is to first check if the swab is at the time of colonisation, or some time later
  // if at the time of colonisation, look for a swab, AND look at the source swabs, and the other source infections
  // if not at the time of colonisation, look for a swab, AND all other infections
  // CONSIDER IMPORTATIONS LATER
  
  IntegerVector t_c = data["t_c"];
  IntegerVector source = data["source"];
  IntegerVector genetic_ids = data["genetic_ids"];
  IntegerVector sample_times = data["sample_times"];
  IntegerVector variant_numbers = data["variant_numbers"];
  IntegerVector variant_idx = WhichVec(current_variant, variant_numbers);
  int col_time = t_c[target];
  
  // First check swabs
  IntegerVector target_genetic_id_locs = WhichVec(target, genetic_ids);
  IntegerVector target_same_variant_locs = intersect(target_genetic_id_locs,variant_idx);
  IntegerVector target_swab_times = as<IntegerVector>(sample_times[target_same_variant_locs]);
  IntegerVector time_loc = WhichVec(time, target_swab_times);
  if(time_loc.length()>0) return true;
  
  // no swab has been found at the time, instead check other individuals

  if(col_time == time)
  {
    // We are checking at the time of colonisation, therefore we wish to first check the source of colonisation
    // and then any other infections at the same time
    int new_target = source[target];
    IntegerVector target_genetic_id_locs = WhichVec(new_target, genetic_ids);
    IntegerVector target_same_variant_locs = intersect(target_genetic_id_locs,variant_idx);
    IntegerVector target_swab_times = as<IntegerVector>(sample_times[target_same_variant_locs]);
    IntegerVector time_loc = WhichVec(time, target_swab_times);
    if(time_loc.length()>0) return true;
    
    // check through other infections
    IntegerVector source_infections = WhichVec(new_target, source);
    for(int i = 0; i<source_infections.length(); i++)
    {
      int current_source_infection = source_infections[i];
      int source_col_time = t_c[current_source_infection];
      if(current_source_infection != target && source_col_time == time)
      {
        new_target = current_source_infection;
        IntegerVector target_genetic_id_locs = WhichVec(new_target, genetic_ids);
        IntegerVector target_same_variant_locs = intersect(target_genetic_id_locs,variant_idx);
        IntegerVector target_swab_times = as<IntegerVector>(sample_times[target_same_variant_locs]);
        IntegerVector time_loc = WhichVec(time, target_swab_times);
        if(time_loc.length()>0) return true;
      }
    }
  }
  else
  {
    // We are checking some time after colonisation, look through other infections at that time
    
    IntegerVector target_infections = WhichVec(target,source);
    for(int i = 0; i<target_infections.length(); i++)
    {
      int current_target = target_infections[i];
      int target_col_time = t_c[current_target];
      if(target_col_time == time)
      {
        // they were colonised at the time we are looking for
        IntegerVector target_genetic_id_locs = WhichVec(current_target, genetic_ids);
        IntegerVector target_same_variant_locs = intersect(target_genetic_id_locs,variant_idx);
        IntegerVector target_swab_times = as<IntegerVector>(sample_times[target_same_variant_locs]);
        IntegerVector time_loc = WhichVec(time, target_swab_times);
        if(time_loc.length()>0) return true;
      }
    }
  }
  return false;
}



/* This function will take the current configuration of the transmission tree, and impute nodes that are needed.
 * This is achieved by sequentially checking each individual if the imputation is needed. If a node is required to be
 * imputed, then the function will check if there is already an imputed node at the current configuration, and if there is a link
 * the distances will be copied. Otherwise distances will be randomly generated. Finally the function will check the previous
 * configuration and find imputed links which are now redundant, and the algorithm will return an updated list and the contribution
 * to the proposal ratio.
 */
List UpdateImputedNodes(List data, NumericVector parameters)
{
  bool print_debug = true;
  double mutation_rate = parameters[3];
  int seq_length = data["sequence_length"];
  // First we wish to use ONLY the observed data, and work out nodes that need to be imputed
  IntegerVector t_c = data["t_c"];
  IntegerVector source = data["source"];
  IntegerVector sample_times = data["sample_times"];
  IntegerVector imputed_nodes = data["imputed_nodes"];
  IntegerVector genetic_ids = data["genetic_ids"];
  IntegerVector variant_numbers = data["variant_numbers"];
  IntegerMatrix genetic_matrix = data["genetic_matrix"];
  
  
  // Load observed data 
  IntegerVector observed_loc = WhichVec(0,imputed_nodes);
  IntegerVector observed_genetic_ids = as<IntegerVector>(genetic_ids[observed_loc]);
  IntegerVector observed_sample_times = as<IntegerVector>(sample_times[observed_loc]);
  IntegerVector observed_variant_numbers = as<IntegerVector>(variant_numbers[observed_loc]);
  IntegerMatrix observed_genetic_matrix(observed_loc.length());
  for(int i = 0; i<observed_loc.length(); i++)
  {
    for(int j = 0; j<observed_loc.length(); j++)
    {
      observed_genetic_matrix(i,j) = genetic_matrix(i,j);
    }
  }
  
  IntegerVector genetic_ids_can = clone(observed_genetic_ids);
  IntegerVector sample_times_can = clone(observed_sample_times);
  IntegerVector variant_numbers_can = clone(observed_variant_numbers);
  IntegerMatrix genetic_matrix_can = clone(observed_genetic_matrix);
  
  IntegerVector gen_source = ReturnGenSourceVector_WHD(data);
  IntegerVector already_seen_links;
  double genetic_contribution = 0;
  
  // Loop through each of the individuals and check if they need a node imputed
  for(int patient_id = 0; patient_id<t_c.length(); patient_id++)
  {
    if(t_c[patient_id] != -1 && source[patient_id] != -1)
    {
      int current_variant = 1; // NEED TO FIX THIS LATER WHEN DEALING WITH MULTIPLE VARIANTS
      int col_time = t_c[patient_id];
      if(print_debug) Rcout << "Check patient " << patient_id << ", with t_c = " << col_time << std::endl;
      IntegerVector patient_swab_loc = WhichVec(patient_id, genetic_ids_can);
      IntegerVector variant_loc = WhichVec(current_variant, variant_numbers_can);
      IntegerVector patient_variant_swab_locs = intersect(patient_swab_loc, variant_loc);
      IntegerVector patient_swab_times = as<IntegerVector>(sample_times_can[patient_variant_swab_locs]);
      IntegerVector patient_swab_loc_same_time = WhichVec(col_time, patient_swab_times);
      
      if(patient_swab_loc_same_time.length()==0)
      {
        if(print_debug) Rcout << "No swab found at the time of infection - checking infector swabs" << std::endl;
        // Now check infector swabs and infector other infections
        int infector_id = source[patient_id];
        IntegerVector infector_swab_loc = WhichVec(infector_id, genetic_ids_can);
        IntegerVector infector_variant_swab_locs = intersect(infector_swab_loc, variant_loc);
        IntegerVector infector_swab_times = as<IntegerVector>(sample_times_can[infector_variant_swab_locs]);
        IntegerVector infector_swab_loc_same_time = WhichVec(col_time, infector_swab_times);
        if(infector_swab_loc_same_time.length()==0)
        {
          if(print_debug) Rcout << "No swab found in the infector, checking their infections for swabs at the time of colonisation" << std::endl;
          // The infector also does not have a swab at the time of infection, check if there are any other infections at this time
          IntegerVector infectees = WhichVec(infector_id, source);
          
          bool swab_found = false;
          for(int j = 0; j<infectees.length(); j++)
          {
            int current_infectee = infectees[j];
            if(current_infectee != patient_id && t_c[current_infectee] == t_c[patient_id])
            {
              // The infectee is not the original target, and also has a colonisation at the same time
              IntegerVector infectee_swab_loc = WhichVec(current_infectee, genetic_ids_can);
              IntegerVector infectee_variant_swab_locs = intersect(infectee_swab_loc, variant_loc);
              IntegerVector infectee_swab_times = as<IntegerVector>(sample_times_can[infectee_variant_swab_locs]);
              IntegerVector infectee_swab_loc_same_time = WhichVec(col_time, infectee_swab_times);
              if(infectee_swab_loc_same_time.length()!=0)
              {
                swab_found = true;
                break;
              }
            }
          }
          
          if(!swab_found)
          {
            if(print_debug) Rcout << "No swab has been found, therefore we must impute distances" << std::endl;
            // there has been no swab found at the same time, therefore try and add a node
            IntegerVector updated_genetic_ids = clone(genetic_ids_can);
            IntegerVector updated_sample_times = clone(sample_times_can);
            IntegerVector updated_variant_numbers = clone(variant_numbers_can);
            updated_genetic_ids.push_back(patient_id);
            updated_sample_times.push_back(col_time);
            updated_variant_numbers.push_back(current_variant);
            if(print_debug)
            {
              Rcout << "Updated genetic ids = " << updated_genetic_ids << std::endl;
              Rcout << "Updated sample times = " << updated_sample_times << std::endl;
              Rcout << "Updated variant numbers = " << updated_variant_numbers << std::endl;
            }
            List temp_data = List::create(Named("t_c") = t_c,
                                          Named("source") = source,
                                          Named("genetic_ids") = updated_genetic_ids,
                                          Named("sample_times") = updated_sample_times,
                                          Named("variant_numbers") = updated_variant_numbers);
            IntegerVector updated_gen_source = ReturnGenSourceVector_WHD(temp_data);
            if(print_debug) Rcout << "Updated genetic source vector = " << updated_gen_source << std::endl;
            IntegerVector children_nodes = WhichVec(updated_gen_source.length()-1, updated_gen_source);
            if(print_debug) Rcout << "Children nodes = " << children_nodes << std::endl;
            if(children_nodes.length() > 1)
            {
              // The imputed node is valid (has 2 or more links, therefore impute matrix)
              int imputed_node = updated_gen_source.length()-1;
              int parent_node = updated_gen_source[imputed_node];
              int imputed_node_cur = -1;
              int target_node = -1;
              
              
              // check if this imputed node has been seen before
              bool already_seen_link = false;
              for(int i = observed_loc.length(); i<genetic_ids.length(); i++)
              {
                if(genetic_ids[i] == patient_id && parent_node == gen_source[i] && sample_times[i] == col_time)
                {
                  already_seen_link = true;
                  imputed_node_cur = i;
                  already_seen_links.push_back(i);
                  break;
                }
              }
              
              // work out the distance between the target and the imputed node
              int draw;
              if(already_seen_link)
              {
                if(print_debug) Rcout << "Already seen the link - copy from previous matrix" << std::endl;
                // the link has already been seen, copy it from the current matrix
                if(parent_node == -1)
                {
                  // the node we are imputing is an imported case, find the one with the smallest time
                  int min_node = children_nodes[0];
                  int min_time = sample_times[min_node];
                  for(int j = 1; j<children_nodes.length(); j++)
                  {
                    int current_child = children_nodes[j];
                    int current_time = sample_times[current_child];
                    if(current_time < min_time)
                    {
                      min_node = current_child;
                      min_time = current_time;
                    }
                  }
                  
                  target_node = min_node;
                  draw = genetic_matrix(imputed_node_cur, min_node);
                }
                else
                {
                  // the node we are imputing is an internal node
                  if(imputed_node_cur == -1) stop("imputed node is -1");
                  draw = genetic_matrix(imputed_node_cur, parent_node);
                  target_node = parent_node;
                }
              }
              else
              {
                // the link has not been seen before, therefore generate the links
                if(print_debug) Rcout << "New link - generate the distances" << std::endl;
                if(parent_node == -1)
                {
                  // We are imputing an imported case
                  // Find the child node with the minimum time
                  int min_node = children_nodes[0];
                  int min_time = updated_sample_times[min_node];
                  for(int j = 1; j<children_nodes.length(); j++)
                  {
                    int current_child = children_nodes[j];
                    int current_time = updated_sample_times[current_child];
                    if(current_time < min_time)
                    {
                      min_node = current_child;
                      min_time = current_time;
                    }
                  }
                  double min_time_double = min_time - updated_sample_times[imputed_node];
                  if(print_debug) Rcout << "min time = " << min_time << ", imputed node time = " << updated_sample_times[imputed_node] << std::endl;
                  double mut_prob = CalculateMutationProbability(mutation_rate, min_time_double); 
                  draw = R::rbinom(seq_length, mut_prob);
                  if(print_debug) Rcout << "Min time double = " << min_time_double << ", mut prob = " << mut_prob << ", draw = " << draw << std::endl;
                  genetic_contribution -= R::dbinom(draw, seq_length, mut_prob, 1);
                  target_node = min_node;
                }
                else
                {
                  // We are imputing an internal node
                  int min_node = children_nodes[0];
                  int min_distance = genetic_matrix(parent_node, children_nodes[0]);
                  for(int j = 1; j<children_nodes.length(); j++)
                  {
                    int current_child = children_nodes[j];
                    int current_distance = genetic_matrix(parent_node, current_child);
                    if(current_distance < min_distance)
                    {
                      min_node = current_child;
                      min_distance = current_distance;
                    }
                  }
                  if(print_debug)
                  {
                    Rcout << "Parent node = " << parent_node << ", min node = " << min_node << ", distance = " << min_distance << std::endl;
                  }
                  // Consider changing to updated_sample_times !
                  int max_time_diff = sample_times_can[min_node] - sample_times_can[parent_node];
                  int current_time_diff = col_time - sample_times_can[parent_node];
                  double time_ratio = (double)(current_time_diff)/(double)(max_time_diff);
                  if(max_time_diff == 0) stop("error in logic, there is a node at the same time");
                  
                  draw = R::rbinom(min_distance, time_ratio);
                  if(draw != 0) Rcout << "Draw = " << draw << std::endl;
                  if(print_debug) Rcout << "Gen contribution before = " << genetic_contribution << std::endl;
                  genetic_contribution -= R::dbinom(draw, min_distance, time_ratio, 1);
                  if(print_debug) Rcout << "Gen contribution after = " << genetic_contribution << std::endl;
                  target_node = parent_node;
                  if(print_debug) Rcout << "draw = " << draw << ", max time diff = " << max_time_diff << ", current time diff = " 
                                        << current_time_diff << ", ratio = " << time_ratio << std::endl;
                }
              }
              
              
              
              // update the genetic matrix
              // initialise a new genetic matrix of zeros
              IntegerMatrix updated_genetic_matrix(updated_gen_source.length());
              
              // copy entries from the old matrix
              for(int ii = 0; ii<updated_gen_source.length()-1; ii++)
              {
                for(int jj = 0; jj<updated_gen_source.length()-1; jj++)
                {
                  updated_genetic_matrix(ii,jj) = genetic_matrix_can(ii,jj);
                }
              }
              
              if(print_debug)
              {
                Rcout << "Printing genetic matrix" << std::endl;
                for(int ii = 0; ii<genetic_matrix.nrow(); ii++)
                {
                  for(int jj = 0; jj<genetic_matrix.ncol(); jj++)
                  {
                    Rcout << genetic_matrix(ii,jj) << " ";
                  }
                  Rcout << std::endl;
                }
              }
              
              if(print_debug)
              {
                Rcout << "Printing genetic matrix can" << std::endl;
                for(int ii = 0; ii<genetic_matrix_can.nrow(); ii++)
                {
                  for(int jj = 0; jj<genetic_matrix_can.ncol(); jj++)
                  {
                    Rcout << genetic_matrix_can(ii,jj) << " ";
                  }
                  Rcout << std::endl;
                }
              }
              
              if(print_debug)
              {
                Rcout << "Printing updated genetic matrix" << std::endl;
                for(int ii = 0; ii<updated_genetic_matrix.nrow(); ii++)
                {
                  for(int jj = 0; jj<updated_genetic_matrix.ncol(); jj++)
                  {
                    Rcout << updated_genetic_matrix(ii,jj) << " ";
                  }
                  Rcout << std::endl;
                }
              }
              
              
              if(print_debug) Rcout << "Imputed node = " << imputed_node << std::endl;
              // now fill in final column
              updated_genetic_matrix(target_node,imputed_node) = draw;
              updated_genetic_matrix(imputed_node,target_node) = draw;
              
              
              // now fill in the rest
              for(int j = 0; j<updated_gen_source.length()-1; j++)
              {
                if(j != target_node)
                {
                  IntegerVector path = ReturnPathToRoot_Rcpp(updated_gen_source, j);
                  if(print_debug) Rcout << "Node = " << j << ", path = " << path << std::endl;
                  bool path_contains_min_node = false;
                  for(int k = 0; k<path.length(); k++)
                  {
                    if(path[k] == target_node)
                    {
                      path_contains_min_node = true;
                      break;
                    }
                  }
                  if(path_contains_min_node)
                  {
                    if(parent_node == -1)
                    {
                      updated_genetic_matrix(j,imputed_node) = updated_genetic_matrix(target_node,j) + updated_genetic_matrix(target_node,imputed_node);
                      updated_genetic_matrix(imputed_node,j) = updated_genetic_matrix(target_node,j) + updated_genetic_matrix(target_node,imputed_node);
                    }
                    else
                    {
                      updated_genetic_matrix(j,imputed_node) = updated_genetic_matrix(target_node,j) - updated_genetic_matrix(target_node,imputed_node);
                      updated_genetic_matrix(imputed_node,j) = updated_genetic_matrix(target_node,j) - updated_genetic_matrix(target_node,imputed_node);
                    }
                    
                    
                  }
                  else
                  {
                    if(parent_node == -1)
                    {
                      updated_genetic_matrix(j,imputed_node) = updated_genetic_matrix(target_node,j) - updated_genetic_matrix(target_node,imputed_node);
                      updated_genetic_matrix(imputed_node,j) = updated_genetic_matrix(target_node,j) - updated_genetic_matrix(target_node,imputed_node);
                    }
                    else
                    {
                      updated_genetic_matrix(j,imputed_node) = updated_genetic_matrix(target_node,j) + updated_genetic_matrix(target_node,imputed_node);
                      updated_genetic_matrix(imputed_node,j) = updated_genetic_matrix(target_node,j) + updated_genetic_matrix(target_node,imputed_node);
                    }
                    
                  }
                }
              }
              if(print_debug)
              {
                Rcout << "Printing updated genetic matrix" << std::endl;
                for(int ii = 0; ii<updated_genetic_matrix.nrow(); ii++)
                {
                  for(int jj = 0; jj<updated_genetic_matrix.ncol(); jj++)
                  {
                    Rcout << updated_genetic_matrix(ii,jj) << " ";
                  }
                  Rcout << std::endl;
                }
              }
              
              genetic_ids_can = updated_genetic_ids;
              sample_times_can = updated_sample_times;
              variant_numbers_can = updated_variant_numbers;
              genetic_matrix_can = updated_genetic_matrix;
              
            }
          }
          
        }
      }
    }
  }
  

  
  int original_length = observed_genetic_ids.length();
  int new_length = genetic_ids_can.length();
  
  // Now calculate previous links that will have to be removed
  for(int i = observed_loc.length(); i<genetic_ids.length(); i++)
  {
    IntegerVector loc_already_seen = WhichVec(i,already_seen_links);
    if(loc_already_seen.length()==0)
    {
      // the link was not reused, therefore calculate the genetic contribution
      int node = i;
      int parent = gen_source[node];
      IntegerVector children_nodes = WhichVec(i, gen_source);
      if(parent == -1)
      {
        // we are removing an imported node
        // the node we are imputing is an imported case, find the one with the smallest time
        int min_node = children_nodes[0];
        int min_time = sample_times[min_node];
        for(int j = 1; j<children_nodes.length(); j++)
        {
          int current_child = children_nodes[j];
          int current_time = sample_times[current_child];
          if(current_time < min_time)
          {
            min_node = current_child;
            min_time = current_time;
          }
        }
        

        int draw = genetic_matrix(i, min_node);
        double min_time_double = min_time - sample_times[node];
        double mut_prob = CalculateMutationProbability(mutation_rate, min_time_double); 
        genetic_contribution += R::dbinom(draw, seq_length, mut_prob, 1);
      }
      else
      {
        // we are removing an internal node
        int parent_node = parent;
        int min_node = children_nodes[0];
        int min_distance = genetic_matrix(parent_node, children_nodes[0]);
        for(int j = 1; j<children_nodes.length(); j++)
        {
          int current_child = children_nodes[j];
          int current_distance = genetic_matrix(parent_node, current_child);
          if(current_distance < min_distance)
          {
            min_node = current_child;
            min_distance = current_distance;
          }
        }
        if(print_debug)
        {
          Rcout << "Parent node = " << parent_node << ", min node = " << min_node << ", distance = " << min_distance << std::endl;
        }
        
        int max_time_diff = sample_times[min_node] - sample_times[parent_node];
        int current_time_diff = sample_times[i] - sample_times[parent_node];
        double time_ratio = (double)(current_time_diff)/(double)(max_time_diff);
        if(max_time_diff == 0) stop("error in logic, there is a node at the same time");
        int draw = genetic_matrix(i,parent_node);
        
        genetic_contribution += R::dbinom(draw, min_distance, time_ratio, 1);
      }
    }
  }
  
  //Rcout << "genetic contribution = " << genetic_contribution << std::endl;
  
  
  IntegerVector imputed_nodes_can(new_length);
  if(print_debug) Rcout << "original length = " << original_length << ", new length = " << new_length << std::endl;
  if(print_debug) Rcout << "Genetic contribution = " << genetic_contribution << std::endl;
  for(int i = original_length; i<new_length; i++)
  {
    imputed_nodes_can[i] = 1;
  }
  
  
  
  if(genetic_ids_can.length() != genetic_ids_can.length() || genetic_ids_can.length() != variant_numbers_can.length()
       || genetic_ids_can.length() != genetic_matrix_can.nrow())
  {
    Rcout << "Gen length = " << genetic_ids_can.length() << ", sample length = " << sample_times_can.length() 
          << ", variant numbers length = " << variant_numbers_can << ", nrow gen mat = " << genetic_matrix_can.nrow() << std::endl;
    stop("inconsistant dimensions for augmented data");
  }
  
  if(print_debug) Rcout << "End of fn - genetic ids = " << genetic_ids << std::endl;
  List return_list = List::create(Named("t_a") = data["t_a"], 
                                  Named("t_c") = t_c,
                                  Named("t_d") = data["t_d"],
                                  Named("source") = source,
                                  Named("screening_matrix") = data["screening_matrix"],
                                  Named("genetic_ids") = genetic_ids_can,
                                  Named("sample_times") = sample_times_can,
                                  Named("variant_numbers") = variant_numbers_can,
                                  Named("genetic_matrix") = genetic_matrix_can,
                                  Named("sequence_length") = data["sequence_length"],
                                  Named("imputed_nodes") = imputed_nodes_can,
                                  Named("genetic_contribution") = genetic_contribution);
  return return_list;
}

bool DoesTargetHaveImputedNodeAtTime(List data, int target, int time)
{
  IntegerVector genetic_ids = data["genetic_ids"];
  IntegerVector imputed_nodes = data["imputed_nodes"];
  IntegerVector sample_times = data["sample_times"];
  IntegerVector imputed_nodes_idx = WhichVec(1, imputed_nodes);
  bool target_has_imputed_node_at_time = false;
  IntegerVector target_genetic_idx = WhichVec(target, genetic_ids);
  IntegerVector target_imputed_genetic_idx = intersect(imputed_nodes_idx, target_genetic_idx);
  

  for(int i = 0; i<target_imputed_genetic_idx.length(); i++)
  {
    int current_node = target_imputed_genetic_idx[i];
    int current_time = sample_times[current_node];
    if(current_time == time)
    {
      target_has_imputed_node_at_time = true;
    }
  }
  
  return target_has_imputed_node_at_time;
}

bool DoesTargetHaveSequenceAtTime(List data, int target, int time)
{
  IntegerVector genetic_ids = data["genetic_ids"];
  IntegerVector sample_times = data["sample_times"];
  IntegerVector genetic_idx = WhichVec(target, genetic_ids);
  // NEED TO INCORPORATE VARIANTS
  IntegerVector target_sample_times = as<IntegerVector>(sample_times[genetic_idx]);
  IntegerVector time_loc = WhichVec(time, target_sample_times);
  if(time_loc.length()==0)
  {
    return false;
  }
  else
  {
    return true;
  }

}


List UpdateImputedNodesMove(List data_cur, List data_can, NumericVector parameters, int target)
{
  IntegerVector imputed_nodes = data_cur["imputed_nodes"];
  IntegerVector imputed_idx = WhichVec(1, imputed_nodes);
  IntegerVector sample_times_cur = data_cur["sample_times"];
  IntegerVector genetic_ids_cur = data_cur["genetic_ids"];
  IntegerVector variant_numbers_cur = data_cur["variant_numbers"];
  IntegerMatrix genetic_matrix_cur = data_cur["genetic_matrix"];
  
  double genetic_contribution = 0;
  
  // First check the candidate gen source and identify any links that don't need to be there
  IntegerVector gen_source_cur = ReturnGenSourceVector_WHD(data_cur);
  IntegerVector gen_source_can = ReturnGenSourceVector_WHD(data_can);
  IntegerVector imputed_nodes_to_remove;
  
  for(int i = 0; i<imputed_idx.length(); i++)
  {
    int current_imputed_node = imputed_idx[i];
    IntegerVector children_nodes = WhichVec(current_imputed_node, gen_source_can);
    if(children_nodes.length()<2)
    {
      // We have found a redundant link, therefore we should calculate the genetic contribution and then remove it
      int parent_node = gen_source_cur[current_imputed_node];
      if(parent_node == -1)
      {
        // We are removing an exterior node
        
        // Look at the current configuration
        IntegerVector children_nodes_cur = WhichVec(current_imputed_node, gen_source_cur);
        // find the minimum of the pairwise distances
        int min_distance = 1000000000;
        int min_i = -1;
        int min_j = -1;
        for(int i = 0; i<(children_nodes_cur.length()-1); i++)
        {
          for(int j = (i+1); j<children_nodes_cur.length(); j++)
          {
            int node_i = children_nodes_cur[i];
            int node_j = children_nodes_cur[j];
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
          Rcout << "Current imputed node = " << current_imputed_node << std::endl;
          Rcout << "Children nodes cur = " << children_nodes_cur << std::endl;
          Rcout << "Gen source cur = " << gen_source_cur << std::endl;
          Rcout << "Genetic ids = " << genetic_ids_cur << std::endl;
          Rcout << "Sample times = " << sample_times_cur << std::endl;
          stop("Problem trying to find the minimum pairwise distance when removing an exterior node");
        }
        
        
        // simulate the distance between min node i and the imputed node
        double upper = sample_times_cur[min_i] - sample_times_cur[current_imputed_node];
        double lower = sample_times_cur[min_i] + sample_times_cur[min_j] - 2*sample_times_cur[current_imputed_node];
        double time_ratio = upper/lower;
        if(time_ratio < 0 || time_ratio > 1)
        {
          Rcout << "Time ratio = " << time_ratio << std::endl;
          stop("invalid time ratio");
        }
        
        int observed_distance = genetic_matrix_cur(current_imputed_node, min_i);
        genetic_contribution -= R::dbinom(observed_distance, min_distance, time_ratio, 1);
        imputed_nodes_to_remove.push_back(current_imputed_node);
      }
      else
      {
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
        double upper = sample_times_cur[current_imputed_node] - sample_times_cur[parent_node];
        double lower = sample_times_cur[min_child] - sample_times_cur[parent_node];
        double time_ratio = (double)(upper)/(double)(lower);
        int observed_distance = genetic_matrix_cur(parent_node, current_imputed_node);
        // NEED TO CHECK THIS - COULD BE MINUS
        genetic_contribution += R::dbinom(observed_distance, min_distance, time_ratio, 1);
        imputed_nodes_to_remove.push_back(current_imputed_node);
      }
    }
  }
  
  // Now update the matrix and vectors
  IntegerVector updated_genetic_ids = clone(genetic_ids_cur);
  IntegerVector updated_sample_times = clone(sample_times_cur);
  IntegerVector updated_imputed_nodes = clone(imputed_nodes);
  IntegerVector updated_variant_numbers = clone(variant_numbers_cur);
  IntegerMatrix updated_genetic_matrix = clone(genetic_matrix_cur);

  
  for(int i = 0; i<imputed_nodes_to_remove.length(); i++)
  {
    int idx_to_remove = imputed_nodes_to_remove[i];
    updated_genetic_ids = RemoveElement(idx_to_remove-i, updated_genetic_ids);
    updated_sample_times = RemoveElement(idx_to_remove-i, updated_sample_times);
    updated_imputed_nodes = RemoveElement(idx_to_remove-i, updated_imputed_nodes);
    updated_variant_numbers = RemoveElement(idx_to_remove-i, updated_variant_numbers);
    updated_genetic_matrix = crossout(updated_genetic_matrix,idx_to_remove-i);
  }
  
  
  if(imputed_nodes_to_remove.length() > 0)
  {
    Rcout << "Nodes " << imputed_nodes_to_remove << " have been removed" << std::endl;
  }
  
  
  // Now impute any new nodes if necessary
  //Rcout << "Attempt to impute nodes" << std::endl;
  IntegerVector source_can = data_can["source"];
  IntegerVector t_c_can = data_can["t_c"];
  // Add a node to the source at the time of colonisation, and check links
  int target_source = source_can[target];
  
  
  // Check if we are proposing an importation, if so exit 
  if(target_source == -1)
  {
    // We are proposing an importation, therefore we dont need to impute any nodes
    IntegerVector t_a = data_can["t_a"];
    IntegerVector t_d = data_can["t_d"];
    IntegerMatrix screening_matrix = data_can["screening_matrix"];
    int sequence_length = data_can["sequence_length"];
    
    
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
                                   Named("genetic_contribution") = genetic_contribution);
    return final_data;
  }
  
  int col_time_can = t_c_can[target];
  //Rcout << "Target = " << target << ", t_c_can = " << col_time_can << ", source = " << target_source << std::endl;
  
  // Check if there is already a sequence in the source at the time
  IntegerVector source_genetic_ids = WhichVec(target_source, updated_genetic_ids);
  IntegerVector target_genetic_ids = WhichVec(target, updated_genetic_ids);
  bool observed_sequence_found = false;
  for(int i = 0; i<source_genetic_ids.length(); i++)
  {
    int current_idx = source_genetic_ids[i];
    int current_time = updated_sample_times[current_idx];
    if(current_time == col_time_can)
    {
      observed_sequence_found = true;
      break;
    }
  }
  
  for(int i = 0; i<target_genetic_ids.length(); i++)
  {
    int current_idx = target_genetic_ids[i];
    int current_time = updated_sample_times[current_idx];
    if(current_time == col_time_can)
    {
      observed_sequence_found = true;
      break;
    }
  }
  
  if(observed_sequence_found)
  {
    // There is already a sequence at that time, dont need to add anything new
    IntegerVector t_a = data_can["t_a"];
    IntegerVector t_d = data_can["t_d"];
    IntegerMatrix screening_matrix = data_can["screening_matrix"];
    int sequence_length = data_can["sequence_length"];
    
    
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
                                   Named("genetic_contribution") = genetic_contribution);
    return final_data;
  }
  
  
  
  IntegerVector genetic_ids_can = clone(updated_genetic_ids);
  IntegerVector sample_times_can = clone(updated_sample_times);
  IntegerVector imputed_nodes_can = clone(updated_imputed_nodes);
  IntegerVector variant_numbers_can = clone(updated_variant_numbers);
  
  genetic_ids_can.push_back(target_source);
  sample_times_can.push_back(col_time_can);
  variant_numbers_can.push_back(1);
  imputed_nodes_can.push_back(1);
  
  List temp_data = List::create(Named("t_c") = t_c_can,
                                Named("source") = source_can,
                                Named("genetic_ids") = genetic_ids_can,
                                Named("sample_times") = sample_times_can,
                                Named("variant_numbers") = variant_numbers_can);
  gen_source_can = ReturnGenSourceVector_WHD(temp_data);
  int new_length = genetic_ids_can.length();
  int imputed_node = new_length - 1;
  IntegerVector children_nodes = WhichVec(imputed_node, gen_source_can);
  IntegerMatrix genetic_matrix_can(new_length);
  if(children_nodes.length()>1)
  {
    // We need to impute a node
    int parent_node = gen_source_can[imputed_node];
    
    // copy the genetic matrix 
    for(int i = 0; i<(new_length-1); i++)
    {
      for(int j = 0; j<(new_length-1); j++)
      {
        genetic_matrix_can(i,j) = updated_genetic_matrix(i,j);
      }
    }
    
    // check if interior or exterior node
    if(parent_node == -1)
    {
      // We need to impute an exterior node
      
      // find the minimum of the pairwise distances
      int min_distance = 1000000000;
      int min_i = -1;
      int min_j = -1;
      Rcout << "Children nodes = " << children_nodes << std::endl;
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
        stop("Problem trying to find the minimum pairwise distance!");
      }
      
      
      // simulate the distance between min node i and the imputed node
      double upper = sample_times_can[min_i] - sample_times_can[imputed_node];
      double lower = sample_times_can[min_i] + sample_times_can[min_j] - 2*sample_times_can[imputed_node];
      double time_ratio = upper/lower;
      if(time_ratio < 0 || time_ratio > 1)
      {
        Rcout << "Time ratio = " << time_ratio << std::endl;
        stop("invalid time ratio");
      }
      
      int simulated_distance = R::rbinom(min_distance, time_ratio);
      genetic_contribution -= R::dbinom(simulated_distance, min_distance, time_ratio, 1);
        
      Rcout << "Min i = " << min_i << ", min j = " << min_j << ", simulated distance = " << simulated_distance << std::endl;  
      genetic_matrix_can(min_i, imputed_node) = simulated_distance;
      genetic_matrix_can(imputed_node, min_i) = simulated_distance;
      for(int i = 0; i<(new_length-1); i++)
      {
        IntegerVector loc_in_child_nodes = WhichVec(i, children_nodes);
        if(i != min_i)
        {
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
      // Impute an interior node
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
      
      double upper = (double)(sample_times_can[imputed_node]) - (double)(sample_times_can[parent_node]);
      double lower = (double)(sample_times_can[min_child]) - (double)(sample_times_can[parent_node]);
      double time_ratio = upper/lower;
      if(time_ratio < 0 || time_ratio > 1)
      {
        Rcout << "Time ratio = " << time_ratio << std::endl;
        stop("invalid time ratio");
      }
      
      // Simulate a distance
      int simulated_distance = R::rbinom(min_distance, time_ratio);
      genetic_contribution -= R::dbinom(simulated_distance, min_distance, time_ratio, 1);
      
      // fill in the rest 
      // ReturnPathToRoot(IntegerVector gen_source, int node)
      for(int i = 0; i<(new_length-1); i++)
      {
        if(i != min_child)
        {
          IntegerVector path_to_root = ReturnPathToRoot_Rcpp(gen_source_can, i);
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
      }
    }
    
    bool negative_matrix = DoesMatrixContainNegativeEntries(genetic_matrix_can);
    if(negative_matrix)
    {
      Rcout << "gen source can = " << gen_source_can << std::endl;
      PrintMatrix(genetic_matrix_can);
      stop("negative values in the matrix!");
    }
    
    IntegerVector t_a = data_can["t_a"];
    IntegerVector t_d = data_can["t_d"];
    IntegerMatrix screening_matrix = data_can["screening_matrix"];
    int sequence_length = data_can["sequence_length"];
    
    //Rcout << "Need to impute" << std::endl;
    List final_data = List::create(Named("t_a") = t_a, 
                                   Named("t_c") = t_c_can,
                                   Named("t_d") = t_d,
                                   Named("source") = source_can,
                                   Named("screening_matrix") = screening_matrix,
                                   Named("genetic_ids") = genetic_ids_can,
                                   Named("sample_times") = sample_times_can,
                                   Named("variant_numbers") = variant_numbers_can,
                                   Named("genetic_matrix") = genetic_matrix_can,
                                   Named("sequence_length") = sequence_length,
                                   Named("imputed_nodes") = imputed_nodes_can,
                                   Named("genetic_contribution") = genetic_contribution);
    
    return final_data;
  }
  else
  {
    // don't need to impute anythign
    
    IntegerVector t_a = data_can["t_a"];
    IntegerVector t_d = data_can["t_d"];
    IntegerMatrix screening_matrix = data_can["screening_matrix"];
    int sequence_length = data_can["sequence_length"];
    
    Rcout << "Don't need to impute" << std::endl;
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
                                   Named("genetic_contribution") = genetic_contribution);
    
    return final_data;
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

List UpdateImputedNodesMove2(List data_cur, List data_can, NumericVector parameters, int target)
{
  bool print_details = false;
  
  IntegerVector imputed_nodes = data_cur["imputed_nodes"];
  IntegerVector imputed_idx = WhichVec(1, imputed_nodes);
  IntegerVector sample_times_cur = data_cur["sample_times"];
  IntegerVector genetic_ids_cur = data_cur["genetic_ids"];
  IntegerVector variant_numbers_cur = data_cur["variant_numbers"];
  IntegerMatrix genetic_matrix_cur = data_cur["genetic_matrix"];
  
  IntegerVector source_can = data_can["source"];
  IntegerVector t_c_can = data_can["t_c"];
  int col_time_can = t_c_can[target];
  double genetic_contribution = 0;
  int target_source = source_can[target];
  
  // First check the candidate gen source and identify any links that don't need to be there
  IntegerVector gen_source_cur = ReturnGenSourceVector_WHD(data_cur);
  
  IntegerVector genetic_ids_can = clone(genetic_ids_cur);
  IntegerVector sample_times_can = clone(sample_times_cur);
  IntegerVector variant_numbers_can = clone(variant_numbers_cur);
  IntegerVector imputed_nodes_can = clone(imputed_nodes);
  
  // mcheck they are not an importation
  if(target_source != -1)
  {
    // Check if there needs to be a node added to the candididate nodes
    IntegerVector source_genetic_ids = WhichVec(target_source, genetic_ids_can);
    IntegerVector target_genetic_ids = WhichVec(target, genetic_ids_can);
    bool observed_sequence_found = false;
    for(int i = 0; i<source_genetic_ids.length(); i++)
    {
      int current_idx = source_genetic_ids[i];
      int current_time = sample_times_can[current_idx];
      if(current_time == col_time_can)
      {
        observed_sequence_found = true;
        break;
      }
    }
    
    for(int i = 0; i<target_genetic_ids.length(); i++)
    {
      int current_idx = target_genetic_ids[i];
      int current_time = sample_times_can[current_idx];
      if(current_time == col_time_can)
      {
        observed_sequence_found = true;
        break;
      }
    }
    
    
    if(!observed_sequence_found)
    {
      genetic_ids_can.push_back(target_source);
      sample_times_can.push_back(col_time_can);
      variant_numbers_can.push_back(1);
      imputed_nodes_can.push_back(1);
    }
  }
  
  

  

  
  List temp_data = List::create(Named("t_c") = t_c_can,
                                Named("source") = source_can,
                                Named("genetic_ids") = genetic_ids_can,
                                Named("sample_times") = sample_times_can,
                                Named("variant_numbers") = variant_numbers_can);
  IntegerVector gen_source_can = ReturnGenSourceVector_WHD(temp_data);
  
  
  //IntegerVector gen_source_can = ReturnGenSourceVector_WHD(data_can);
  IntegerVector imputed_nodes_to_remove;
  
  for(int i = 0; i<imputed_idx.length(); i++)
  {
    int current_imputed_node = imputed_idx[i];
    IntegerVector children_nodes = WhichVec(current_imputed_node, gen_source_can);
    if(children_nodes.length()<2)
    {
      // We have found a redundant link, therefore we should calculate the genetic contribution and then remove it
      int parent_node = gen_source_cur[current_imputed_node];
      if(parent_node == -1)
      {
        // We are removing an exterior node
        
        // Look at the current configuration
        IntegerVector children_nodes_cur = WhichVec(current_imputed_node, gen_source_cur);
        // find the minimum of the pairwise distances
        int min_distance = 1000000000;
        int min_i = -1;
        int min_j = -1;
        for(int i = 0; i<(children_nodes_cur.length()-1); i++)
        {
          for(int j = (i+1); j<children_nodes_cur.length(); j++)
          {
            int node_i = children_nodes_cur[i];
            int node_j = children_nodes_cur[j];
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
          Rcout << "Current imputed node = " << current_imputed_node << std::endl;
          Rcout << "Children nodes cur = " << children_nodes_cur << std::endl;
          Rcout << "Gen source cur = " << gen_source_cur << std::endl;
          Rcout << "Gen source can = " << gen_source_can << std::endl;
          Rcout << "Genetic ids = " << genetic_ids_cur << std::endl;
          Rcout << "Sample times = " << sample_times_cur << std::endl;
          stop("Problem trying to find the minimum pairwise distance when removing an exterior node");
        }
        
        
        // simulate the distance between min node i and the imputed node
        double upper = sample_times_cur[min_i] - sample_times_cur[current_imputed_node];
        double lower = sample_times_cur[min_i] + sample_times_cur[min_j] - 2*sample_times_cur[current_imputed_node];
        double time_ratio = upper/lower;
        if(time_ratio < 0 || time_ratio > 1)
        {
          Rcout << "Time ratio = " << time_ratio << std::endl;
          stop("invalid time ratio");
        }
        
        int observed_distance = genetic_matrix_cur(current_imputed_node, min_i);
        genetic_contribution -= R::dbinom(observed_distance, min_distance, time_ratio, 1);
        imputed_nodes_to_remove.push_back(current_imputed_node);
      }
      else
      {
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
        double upper = sample_times_cur[current_imputed_node] - sample_times_cur[parent_node];
        double lower = sample_times_cur[min_child] - sample_times_cur[parent_node];
        double time_ratio = (double)(upper)/(double)(lower);
        int observed_distance = genetic_matrix_cur(parent_node, current_imputed_node);
        // NEED TO CHECK THIS - COULD BE MINUS
        genetic_contribution += R::dbinom(observed_distance, min_distance, time_ratio, 1);
        imputed_nodes_to_remove.push_back(current_imputed_node);
      }
    }
  }
  
  if(print_details) Rcout << "Nodes to remove = " << imputed_nodes_to_remove << std::endl;
  
  // Now update the matrix and vectors
  IntegerVector updated_genetic_ids = clone(genetic_ids_cur);
  IntegerVector updated_sample_times = clone(sample_times_cur);
  IntegerVector updated_imputed_nodes = clone(imputed_nodes);
  IntegerVector updated_variant_numbers = clone(variant_numbers_cur);
  IntegerMatrix updated_genetic_matrix = clone(genetic_matrix_cur);
  
  
  for(int i = 0; i<imputed_nodes_to_remove.length(); i++)
  {
    int idx_to_remove = imputed_nodes_to_remove[i];
    updated_genetic_ids = RemoveElement(idx_to_remove-i, updated_genetic_ids);
    updated_sample_times = RemoveElement(idx_to_remove-i, updated_sample_times);
    updated_imputed_nodes = RemoveElement(idx_to_remove-i, updated_imputed_nodes);
    updated_variant_numbers = RemoveElement(idx_to_remove-i, updated_variant_numbers);
    updated_genetic_matrix = crossout(updated_genetic_matrix,idx_to_remove-i);
  }
  
  
  if(imputed_nodes_to_remove.length() > 0)
  {
    //Rcout << "Nodes " << imputed_nodes_to_remove << " have been removed" << std::endl;
  }
  
  
  // Now impute any new nodes if necessary
  //Rcout << "Attempt to impute nodes" << std::endl;
  //IntegerVector source_can = data_can["source"];
  //IntegerVector t_c_can = data_can["t_c"];
  // Add a node to the source at the time of colonisation, and check links
  
  
  
  // Check if we are proposing an importation, if so exit 
  if(target_source == -1)
  {
    // We are proposing an importation, therefore we dont need to impute any nodes
    IntegerVector t_a = data_can["t_a"];
    IntegerVector t_d = data_can["t_d"];
    IntegerMatrix screening_matrix = data_can["screening_matrix"];
    int sequence_length = data_can["sequence_length"];
    
    
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
                                   Named("genetic_contribution") = genetic_contribution);
    return final_data;
  }
  
  
  //Rcout << "Target = " << target << ", t_c_can = " << col_time_can << ", source = " << target_source << std::endl;
  
  // Check if there is already a sequence in the source at the time
  IntegerVector source_genetic_ids = WhichVec(target_source, updated_genetic_ids);
  IntegerVector target_genetic_ids = WhichVec(target, updated_genetic_ids);
  bool observed_sequence_found = false;
  for(int i = 0; i<source_genetic_ids.length(); i++)
  {
    int current_idx = source_genetic_ids[i];
    int current_time = updated_sample_times[current_idx];
    if(current_time == col_time_can)
    {
      observed_sequence_found = true;
      break;
    }
  }
  
  for(int i = 0; i<target_genetic_ids.length(); i++)
  {
    int current_idx = target_genetic_ids[i];
    int current_time = updated_sample_times[current_idx];
    if(current_time == col_time_can)
    {
      observed_sequence_found = true;
      break;
    }
  }
  
  if(observed_sequence_found)
  {
    // There is already a sequence at that time, dont need to add anything new
    IntegerVector t_a = data_can["t_a"];
    IntegerVector t_d = data_can["t_d"];
    IntegerMatrix screening_matrix = data_can["screening_matrix"];
    int sequence_length = data_can["sequence_length"];
    
    
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
                                   Named("genetic_contribution") = genetic_contribution);
    return final_data;
  }
  
  
  
  genetic_ids_can = clone(updated_genetic_ids);
  sample_times_can = clone(updated_sample_times);
  imputed_nodes_can = clone(updated_imputed_nodes);
  variant_numbers_can = clone(updated_variant_numbers);
  
  genetic_ids_can.push_back(target_source);
  sample_times_can.push_back(col_time_can);
  variant_numbers_can.push_back(1);
  imputed_nodes_can.push_back(1);
  
  temp_data = List::create(Named("t_c") = t_c_can,
                                Named("source") = source_can,
                                Named("genetic_ids") = genetic_ids_can,
                                Named("sample_times") = sample_times_can,
                                Named("variant_numbers") = variant_numbers_can);
  gen_source_can = ReturnGenSourceVector_WHD(temp_data);
  int new_length = genetic_ids_can.length();
  int imputed_node = new_length - 1;
  IntegerVector children_nodes = WhichVec(imputed_node, gen_source_can);
  IntegerMatrix genetic_matrix_can(new_length);
  if(children_nodes.length()>1)
  {
    // We need to impute a node
    int parent_node = gen_source_can[imputed_node];
    
    // copy the genetic matrix 
    for(int i = 0; i<(new_length-1); i++)
    {
      for(int j = 0; j<(new_length-1); j++)
      {
        genetic_matrix_can(i,j) = updated_genetic_matrix(i,j);
      }
    }
    
    // check if interior or exterior node
    if(parent_node == -1)
    {
      if(print_details) Rcout << "Try to impute exterior node = " << std::endl;
      // We need to impute an exterior node
      
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
        stop("Problem trying to find the minimum pairwise distance!");
      }
      
      
      // simulate the distance between min node i and the imputed node
      double upper = sample_times_can[min_i] - sample_times_can[imputed_node];
      double lower = sample_times_can[min_i] + sample_times_can[min_j] - 2*sample_times_can[imputed_node];
      double time_ratio = upper/lower;
      if(time_ratio < 0 || time_ratio > 1)
      {
        Rcout << "Time ratio = " << time_ratio << std::endl;
        stop("invalid time ratio");
      }
      
      int simulated_distance = R::rbinom(min_distance, time_ratio);
      genetic_contribution -= R::dbinom(simulated_distance, min_distance, time_ratio, 1);
      
      //Rcout << "Min i = " << min_i << ", min j = " << min_j << ", simulated distance = " << simulated_distance << std::endl;  
      genetic_matrix_can(min_i, imputed_node) = simulated_distance;
      genetic_matrix_can(imputed_node, min_i) = simulated_distance;
      for(int i = 0; i<(new_length-1); i++)
      {
        IntegerVector loc_in_child_nodes = WhichVec(i, children_nodes);
        if(i != min_i)
        {
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
      // Impute an interior node
      if(print_details) Rcout << "Try to impute interior node" << std::endl;
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
      
      double upper = (double)(sample_times_can[imputed_node]) - (double)(sample_times_can[parent_node]);
      double lower = (double)(sample_times_can[min_child]) - (double)(sample_times_can[parent_node]);
      double time_ratio = upper/lower;
      if(print_details)
      {
        Rcout << "Genetic ids can = " << genetic_ids_can << std::endl;
        Rcout << "Sample times can = " << sample_times_can << std::endl;
        Rcout << "Gen source can = " << gen_source_can << std::endl;
        Rcout << "Parent = " << parent_node << ", imputed node = " << imputed_node << ", min child = " << min_child << "col time can = " << col_time_can << std::endl;
        Rcout << "Imputed node sample time = " << sample_times_can[imputed_node] << ", parent sampel time = " << sample_times_can[parent_node]
              << ", min child sample time = " << sample_times_can[min_child] << std::endl;
      }
      if(time_ratio < 0 || time_ratio > 1)
      {
        Rcout << "Time ratio = " << time_ratio << std::endl;
        stop("invalid time ratio");
      }
      
      // Simulate a distance
      int simulated_distance = R::rbinom(min_distance, time_ratio);
      genetic_contribution -= R::dbinom(simulated_distance, min_distance, time_ratio, 1);
      if(print_details) Rcout << "Simulated distance = " << simulated_distance << ", min distance = " << min_distance
                              << ", time ratio = " << time_ratio << ", genetic contribution " << genetic_contribution << std::endl;
      
      
      // fill in the rest 
      // ReturnPathToRoot(IntegerVector gen_source, int node)
      for(int i = 0; i<(new_length-1); i++)
      {
        if(i != min_child)
        {
          IntegerVector path_to_root = ReturnPathToRoot_Rcpp(gen_source_can, i);
          IntegerVector imputed_loc_in_path = WhichVec(imputed_node, path_to_root);
          if(i == 4)
          {
            // debugging
            //Rcout << "Path to root = " << path_to_root << std::endl;
            //Rcout << "Min child = " << min_child << ", min distance = " << min_distance << ", simulated distance = " << simulated_distance << std::endl;
          }
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
      }
    }
    
    bool negative_matrix = DoesMatrixContainNegativeEntries(genetic_matrix_can);
    if(negative_matrix)
    {
      //Rcout << "gen source can = " << gen_source_can << std::endl;
      //PrintMatrix(genetic_matrix_can);
      
      //Rcout << "Children nodes = " << children_nodes << std::endl;
      Rcout << "Negative entries in the distance matrix - no move is made!" << std::endl;
      genetic_contribution -= 1000000000;
      //stop("negative values in the matrix!");
    }
    
    IntegerVector t_a = data_can["t_a"];
    IntegerVector t_d = data_can["t_d"];
    IntegerMatrix screening_matrix = data_can["screening_matrix"];
    int sequence_length = data_can["sequence_length"];
    
    //Rcout << "Need to impute" << std::endl;
    List final_data = List::create(Named("t_a") = t_a, 
                                   Named("t_c") = t_c_can,
                                   Named("t_d") = t_d,
                                   Named("source") = source_can,
                                   Named("screening_matrix") = screening_matrix,
                                   Named("genetic_ids") = genetic_ids_can,
                                   Named("sample_times") = sample_times_can,
                                   Named("variant_numbers") = variant_numbers_can,
                                   Named("genetic_matrix") = genetic_matrix_can,
                                   Named("sequence_length") = sequence_length,
                                   Named("imputed_nodes") = imputed_nodes_can,
                                   Named("genetic_contribution") = genetic_contribution);
    
    return final_data;
  }
  else
  {
    // don't need to impute anythign
    
    IntegerVector t_a = data_can["t_a"];
    IntegerVector t_d = data_can["t_d"];
    IntegerMatrix screening_matrix = data_can["screening_matrix"];
    int sequence_length = data_can["sequence_length"];
    
    //Rcout << "Don't need to impute" << std::endl;
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
                                   Named("genetic_contribution") = genetic_contribution);
    
    return final_data;
  }
  
  
} 
// Could be a useful function
//List RemoveUnnecessaryImputedNodes(List data, List data_cur)



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
  




 
/* Consider imputing nodes and the corresponding forwards and reverse moves
 * For the reverse move, check if they (their source) had an imputed node at their current time of colonisation, if so check if its still needed
 * where we check if there are any other infections at that time, and if it is even needed at all
 */

List UpdateImputedNodesMove3(List data_cur, List data_can, NumericVector parameters, int target)
{
  IntegerVector imputed_nodes_cur = data_cur["imputed_nodes"];
  IntegerVector genetic_ids = data_cur["genetic_ids"];
  IntegerVector sample_times = data_cur["sample_times"];
  Rcout << "Imputed nodes cur = " << imputed_nodes_cur << std::endl;
  
  //IntegerVector nodes_ne
  List nodes_to_impute_cur = ReturnNodesToImpute(data_cur);
  List nodes_to_impute_can = ReturnNodesToImpute(data_can);
  
  IntegerVector nodes_cur = nodes_to_impute_cur["nodes"];
  IntegerVector times_cur = nodes_to_impute_cur["times"];
  IntegerVector nodes_can = nodes_to_impute_can["nodes"];
  IntegerVector times_can = nodes_to_impute_can["times"];
  
  Rcout << "Nodes cur = " << nodes_cur << std::endl;
  Rcout << "Times cur = " << times_cur << std::endl;
  Rcout << "Nodes can = " << nodes_can << std::endl;
  Rcout << "Times can = " << times_can << std::endl;
  
  IntegerVector gen_source_cur = ReturnGenSourceVector_WHD(data_cur);
  Rcout << "Gen source cur = " << gen_source_cur << std::endl;
  
  IntegerVector t_c_can = data_can["t_c"];
  IntegerVector source_can = data_can["source"];
  IntegerVector variant_numbers = data_cur["variant_numbers"];
  
  IntegerVector observed_idx = WhichVec(0, imputed_nodes_cur);
  IntegerVector genetic_ids_can = as<IntegerVector>(genetic_ids[observed_idx]);
  IntegerVector sample_times_can = as<IntegerVector>(sample_times[observed_idx]);
  IntegerVector variant_numbers_can = as<IntegerVector>(variant_numbers[observed_idx]);
  
  
  double genetic_contribution;
  for(int i = 0; i<nodes_can.length(); i++)
  {
    genetic_ids_can.push_back(nodes_can[i]);
    sample_times_can.push_back(times_can[i]);
    variant_numbers_can.push_back(1);
  }
  Rcout << "genetic_ids_can = " << genetic_ids_can << std::endl;
  Rcout << "sample_times_can = " << sample_times_can << std::endl;
  Rcout << "variant_numbers_can = " << variant_numbers_can << std::endl;
  List temp_data = List::create(Named("t_c") = t_c_can,
                                Named("source") = source_can,
                                Named("genetic_ids") = genetic_ids_can,
                                Named("sample_times") = sample_times_can,
                                Named("variant_numbers") = variant_numbers_can);
  IntegerVector gen_source_can = ReturnGenSourceVector_WHD(temp_data);
  
  
  Rcout << "Gen source can = " << gen_source_can << std::endl;
  

  
  // identify those that are (i) the same, (ii) need to be removed and (iii) need to be added
  // then fix the master distances and update those
  
  // first identify links that need to be removed and calculate the genetic contribution
  // then add the links - if the link is previous been seen then keep that distance and no need to update the contribution
  // if the link hasn't been previously seen, then simualte the distance
  
  IntegerVector genetic_ids_to_remove;
  IntegerVector times_of_genetic_ids_to_remove;
  
  for(int i = 0; nodes_cur.length(); i++)
  {
    int current_node = nodes_cur[i];
    int current_time = times_cur[i];
    IntegerVector node_loc_in_can = WhichVec(current_node, nodes_can);
    bool node_in_can_found = false;
    for(int j = 0; j<node_loc_in_can.length(); j++)
    {
      int current_node_loc_in_can = node_loc_in_can[j];
      int time_in_can = times_can[current_node_loc_in_can];
      if(current_time == time_in_can)
      {
        node_in_can_found = true;
      }
    }
    if(!node_in_can_found)
    {
      genetic_ids_to_remove.push_back(current_node);
      times_of_genetic_ids_to_remove.push_back(current_node);
    }
  }
  
  IntegerMatrix genetic_matrix_cur = data_cur["genetic_matrix"];
  // now calculate the contribution of the genetic id removals
  for(int i = 0; i<genetic_ids_to_remove.length(); i++)
  {
    int current_id = genetic_ids_to_remove[i];
    int current_time = times_of_genetic_ids_to_remove[i];
    IntegerVector id_idx = WhichVec(current_id, genetic_ids);
    IntegerVector time_idx = WhichVec(current_time, sample_times);
    IntegerVector both_idx = intersect(id_idx, time_idx);
    if(both_idx.length() != 1)
    {
      stop("duplicate genetic ids when imputing");
    }
    int imputed_node = both_idx[0];
    int parent_node = gen_source_cur[imputed_node];
    IntegerVector children_nodes = WhichVec(imputed_node, gen_source_cur);
    if(parent_node == -1)
    {
      // We are removing an external node, determine the genetic contribution
      int min_distance = 100000000;
      
      int min_i = -1;
      int min_j = -1;
      //Rcout << "Children nodes = " << children_nodes << std::endl;
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
        stop("Problem trying to find the minimum pairwise distance!");
      }
      
      
      // simulate the distance between min node i and the imputed node
      double upper = sample_times_can[min_i] - sample_times_can[imputed_node];
      double lower = sample_times_can[min_i] + sample_times_can[min_j] - 2*sample_times_can[imputed_node];
      double time_ratio = upper/lower;
      if(time_ratio < 0 || time_ratio > 1)
      {
        Rcout << "Time ratio = " << time_ratio << std::endl;
        stop("invalid time ratio");
      }
      
      int observed_distance = genetic_matrix_cur(imputed_node, min_i);
      genetic_contribution += R::dbinom(observed_distance, min_distance, time_ratio, 1);
    }
    else
    {
      // we are removing an internal node
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
      double upper = sample_times[imputed_node] - sample_times[parent_node];
      double lower = sample_times[min_child] - sample_times[parent_node];
      double time_ratio = (double)(upper)/(double)(lower);
      int observed_distance = genetic_matrix_cur(parent_node, imputed_node);
      // NEED TO CHECK THIS - COULD BE MINUS
      genetic_contribution += R::dbinom(observed_distance, min_distance, time_ratio, 1);
    }
  }
  
  
  // now build up the matrix again
  // go through the links, if seen previously then just copy, otherwise simulate
  IntegerVector updated_genetic_ids = as<IntegerVector>(genetic_ids[observed_idx]);
  IntegerVector updated_sample_times = as<IntegerVector>(sample_times[observed_idx]);
  IntegerVector updated_variant_numbers = as<IntegerVector>(variant_numbers[observed_idx]);
  for(int i = 0; i<nodes_can.length(); i++)
  {
    int current_id = nodes_can[i];
    int current_time = times_can[i];
    updated_genetic_ids.push_back(current_id);
    updated_sample_times.push_back(current_time);
    updated_variant_numbers.push_back(1);
    List temp_data = List::create(Named("t_c") = t_c_can,
                                  Named("source") = source_can,
                                  Named("genetic_ids") = updated_genetic_ids,
                                  Named("sample_times") = updated_sample_times,
                                  Named("variant_numbers") = updated_variant_numbers);
    IntegerVector gen_source = ReturnGenSourceVector_WHD(temp_data);
    
    int imputed_node = gen_source.length()-1;
    int parent_node = gen_source[imputed_node];
    
    // Check if the id and time has previously been seen with the same parent
    IntegerVector id_idx = WhichVec(current_id, genetic_ids);
    IntegerVector time_idx = WhichVec(current_time, sample_times);
    IntegerVector both_idx = intersect(id_idx, time_idx);
    
    
    if(parent_node == -1)
    {
      // we are imputing an exterior node
    }
    else
    {
      // we are imputing an interior node
      if(both_idx.length()==0)
      {
        // there is not a sequence, need to simulate distances
      }
      else
      {
        
      }
    }
    
    
    
  }
  
  
  
  return data_cur;
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
  bool print_debug_file = true;
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
    
    
    IntegerVector gen_source = ReturnGenSourceVector_WHD(temp_data);
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
    IntegerVector gen_source = ReturnGenSourceVector_WHD(temp_data);
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
  
  IntegerVector gen_source_cur = ReturnGenSourceVector_WHD(data_cur);
  IntegerVector gen_source_can = ReturnGenSourceVector_WHD(temp_data);
  
  if(print_details) Rcout << "gen source cur = " << gen_source_cur << std::endl;
  if(print_details) Rcout << "gen source can = " << gen_source_can << std::endl;
  
  
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
                                 Named("genetic_contribution") = genetic_contribution);
  
  return final_data;
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

List CalculateImportDistance(List data, NumericVector parameters) {
  IntegerVector gen_source = ReturnGenSourceVector_WHD(data);
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



// [[Rcpp::export]]
void MCMC(List MCMC_options, int sequence_length, IntegerVector t_a, IntegerVector t_c, IntegerVector t_d, IntegerVector source, IntegerMatrix screening_matrix,
          IntegerVector genetic_ids, IntegerVector sample_times, IntegerVector variant_numbers, IntegerMatrix genetic_matrix)
{
  // Initialise data
  List data = InitialiseData(sequence_length, t_a, t_c, t_d, source, screening_matrix, genetic_ids, sample_times, variant_numbers, genetic_matrix);
  
  
  // Initialise Imputed nodes
  data = InitialiseMCMCImputedNodes2(data);
  

    
  
  // Load MCMC options
  const int iterations = MCMC_options["iterations"];
  const int population_size = t_a.length();
  const int num_updates = MCMC_options["num_updates"];
  bool print_matrix = MCMC_options["print_matrix"];
  NumericVector prior_parameters = MCMC_options["prior_parameters"];
  NumericVector initial_chain_state = MCMC_options["initial_chain_state"];
  NumericVector parameters = clone(initial_chain_state);
  NumericVector proposal_variance = MCMC_options["proposal_variance"];
  IntegerVector debug_flags = MCMC_options["debug_flags"];
  std::string output_file = MCMC_options["output_file"];
  std::string coltime_file = MCMC_options["coltime_file"];
  std::string source_file = MCMC_options["source_file"];
  std::string loglik_file = MCMC_options["loglik_file"];
  std::string matrix_file = MCMC_options["matrix_file"];
  std::string augmented_data_file = "augmented_data.dat";
  
  // Impute any missing nodes
  //IntegerVector genetic_ids2 = data["genetic_ids"];
  //Rcout << "Genetic ids = " << genetic_ids2 << std::endl;
  
  //data = ImputeNodes(data, parameters);
  //data = ImputeNodes4(data, parameters);
  //data = UpdateImputedNodes(data, parameters);
  
  bool print_genetic_information = false;
  
  if(print_genetic_information)
  {
    IntegerVector genetic_ids3 = data["genetic_ids"];
    Rcout << "Genetic ids = " << genetic_ids3 << std::endl;
    
    IntegerVector sample_times3 = data["sample_times"];
    Rcout << "sample times = " << sample_times3 << std::endl;
    
    IntegerVector imputed_nodes = data["imputed_nodes"];
    Rcout << "imputed nodes = " << imputed_nodes << std::endl;
    
    IntegerVector gen_source = ReturnGenSourceVector_WHD(data);
    Rcout << "Gen source = " << gen_source << std::endl;
    
    IntegerMatrix gen_mat = data["genetic_matrix"];
    PrintMatrix(gen_mat);
  }

  //return;
  //for(int i = 0; i<gen_mat.nrow(); i++)
  //{
  //  for(int j = 0; j<gen_mat.ncol(); j++)
  //  {
  //    if(gen_mat(i,j)<0) stop("negative values!");
  //  }
  //}
  
  //return;
  
  
  // Calculate other variables
  int true_positives = CalculateTruePositives(data);
  int nacc_beta = 0;
  int nacc_add = 0;
  int nacc_move = 0;
  int nacc_remove = 0;
  //int nacc_source_change = 0;
  //int nacc_swap = 0;
  int nacc_lambda = 0;
  //int nacc_import_swap = 0;
  double loglik = CalculateLogLikelihood(data, parameters);
  Rcout << "Log likelihood calculcated" << std::endl;
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
  
  // Colonisation time sum file
  remove(coltime_file.c_str());
  std::ofstream myfile2;
  myfile2.open(coltime_file.c_str());
  assert(myfile2.is_open());
  PrintColonisationTimeSumToFile(myfile2, data);
  
  // Source file
  remove(source_file.c_str());
  std::ofstream myfile3;
  myfile3.open(source_file.c_str());
  assert(myfile3.is_open());
  PrintIntVectorToFile(myfile3, data["source"]);
  
  // loglik file
  remove(loglik_file.c_str());
  std::ofstream myfile4;
  myfile4.open(loglik_file.c_str());
  assert(myfile4.is_open());
  myfile4 << loglik << std::endl;
  
  // matrix file
  remove(matrix_file.c_str());
  std::ofstream myfile5;
  myfile5.open(matrix_file.c_str());
  assert(myfile5.is_open());
  
  if(print_matrix)
  {
    IntegerMatrix genetic_matrix_cur = data["genetic_matrix"];
    PrintMatrixToFile(myfile5, genetic_matrix_cur);
  }
  
  
  
  
  // augmented data file
  remove(augmented_data_file.c_str());
  std::ofstream myfile6;
  myfile6.open(augmented_data_file.c_str());
  assert(myfile6.is_open());

  bool print_aug_data = false;
  if(print_aug_data)
  {
    PrintIntVectorToFile(myfile6, data["t_c"]);
    PrintIntVectorToFile(myfile6, data["source"]);
    PrintIntVectorToFile(myfile6, data["genetic_ids"]);
    PrintIntVectorToFile(myfile6, data["sample_times"]);
    PrintIntVectorToFile(myfile6, data["variant_numbers"]);
    myfile6 << std::endl;
  }



  
  
  
  
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
      int import_sum = CalculateImportSum(data);
      parameters[1] = R::rbeta(prior_parameters[2] + import_sum, prior_parameters[3] + population_size - import_sum);
    }
    
    loglik = CalculateLogLikelihood(data, parameters);
    
    // Update beta by metropolis hastings
    if(debug_flags[2]==0)
    {
      UpdateTransmissionRate(data, parameters, loglik, proposal_variance[0], prior_parameters[4], nacc_beta);
    }
    
    if(debug_flags[3]==0)
    {
      UpdateMutationRate(data, parameters, loglik, proposal_variance[1], prior_parameters[5], nacc_lambda);
    }
    
    
    List import_distance_info = CalculateImportDistance(data, parameters);
    double distance_sum = import_distance_info["distance_sum"];
    double import_counter = import_distance_info["import_counter"];
    //double variant_distance_sum = import_distance_info["variant_distance_sum"];
    //double variant_import_counter = import_distance_info["variant_import_counter"];
    
    //double gamma_draw = R::rgamma(prior_parameters[6] + distance_sum, 1/(prior_parameters[7] + import_counter));
    //Rcout << "Distance sum = " << distance_sum << ", import counter = " << import_counter << ", gamma draw = " << std::endl;
    //Rcout << "Shape = " << prior_parameters[6] + distance_sum << ", scale = " << 1/(prior_parameters[7] + import_counter) << std::endl;
    //parameters[4] = gamma_draw;
    NumericVector gamma_draws = Rcpp::rgamma(1, prior_parameters[6] + distance_sum, 1/(prior_parameters[7] + import_counter));
    parameters[4] = gamma_draws[0];
    
    /*
    Rcout << "Distance sum = " << distance_sum << ", import counter = " << import_counter << ", gamma draw = " << std::endl;
    Rcout << "Shape = " << prior_parameters[6] + distance_sum << ", scale = " << 1/(prior_parameters[7] + import_counter) << std::endl;
    Rcout << "Draw = " << gamma_draws[0] << std::endl;
    */
    //Rcout << "Rcpp r gamma " << std::endl;
    
    //parameters[5] = R::rgamma(prior_parameters[6] + variant_distance_sum, 1/(prior_parameters[7] + variant_import_counter));
    
    
    
    // Augmented data updates
    if(debug_flags[4]==0)
    {
      for(int update_counter = 0; update_counter < num_updates; update_counter++)
      {
        double w = 0.3;
        double move = R::runif(0.0,3.0);
        augmented_moves_proposed[floor(move)]++;
        if(floor(move) < 1)
        {
          // Move a colonisation time
          //MoveColonisationTime(data, parameters, loglik, w, Vq, Va, nacc_move); 
          
          // try moving col time and genetic data
          data = MoveColonisationTimeGenetic(data, parameters, loglik, w, Vq, Va, nacc_move);

          //Rcout << "Before move" << std::endl;
          //data = MoveColonisationTime_List(data, parameters, loglik, w, Vq, Va, nacc_move);
          //Rcout << "After move" << std::endl;
          //return;
          bool print_details = false;
          if(print_details)
          {
            IntegerVector genetic_ids = data["genetic_ids"];
            IntegerVector sample_times = data["sample_times"];
            IntegerVector imputed_nodes = data["imputed_nodes"];
            IntegerVector gen_source_vec = ReturnGenSourceVector_WHD(data);
            Rcout << "Genetic ids = " << genetic_ids << std::endl;
            Rcout << "Sample times = " << sample_times << std::endl;
            Rcout << "imputed nodes = " << imputed_nodes << std::endl;
            Rcout << "Gen source = " << gen_source_vec << std::endl;
          }
          
          //if(update_counter == 15) return;
          
          if(print_matrix)
          {
            IntegerMatrix genetic_matrix_cur = data["genetic_matrix"];
            PrintMatrixToFile(myfile5, genetic_matrix_cur);
          }
        }
        else if(floor(move) < 2)
        {
          // Add a colonsation time
          //Rcout << "Before add" << std::endl;
          AddColonisationTime(data, parameters, loglik, w, Vs, Va, nacc_add);
          //data = AddColonisationTime_List(data, parameters, loglik, w, Vs, Va, nacc_add);
          //Rcout << "After add" << std::endl;
        }
        else if(floor(move) < 3)
        {
          // Remove a colonisation time
          //Rcout << "Before remove" << std::endl;
          RemoveColonisationTime(data, parameters, loglik, w, Vs, Va, nacc_remove);
          //data = RemoveColonisationTime_List(data, parameters, loglik, w, Vs, Va, nacc_remove);
          //Rcout << "After remove" << std::endl;
        }
        
        
        if(print_aug_data)
        {
          PrintIntVectorToFile(myfile6, data["t_c"]);
          PrintIntVectorToFile(myfile6, data["source"]);
          PrintIntVectorToFile(myfile6, data["genetic_ids"]);
          PrintIntVectorToFile(myfile6, data["sample_times"]);
          PrintIntVectorToFile(myfile6, data["variant_numbers"]);
          myfile6 << std::endl;
        }
        //Rcout << "Update " << update_counter << " completed" << std::endl;

      }
      // Write colonisation time sum to file
      PrintColonisationTimeSumToFile(myfile2, data);
      PrintIntVectorToFile(myfile3, data["source"]);
      
      
    }
    
    //stop("stop early");
    
    
    
    // Write parameters to file
    PrintNumVectorToFile(myfile, parameters);
    myfile4 << loglik << std::endl;
    
    //Rcout << "Iteration " << i << " completed" << std::endl;
    if(i%100==0)
    {
      Rcout << i << std::endl;
    }
    
  }
  
  // Close files
  myfile.close();
  myfile2.close();
  myfile3.close();
  myfile4.close();
  myfile5.close();
  myfile6.close();
  

  double beta_prob = (double)nacc_beta/(double)iterations;
  double lambda_prob = (double)nacc_lambda/(double)iterations;
  double move_prob = (double)nacc_move/(double)(augmented_moves_proposed[0]);
  double add_prob = (double)nacc_add/(double)(augmented_moves_proposed[1]);
  double remove_prob = (double)nacc_remove/(double)(augmented_moves_proposed[2]);
  Rcout << "Beta acceptance = " << beta_prob << ", Lambda acceptance = " << lambda_prob << std::endl;
  Rcout << "Move prob = " << move_prob << ", add prob = " << add_prob << ", remove prob = " << remove_prob << std::endl;
  
  /*
  double source_change_prob = (double)nacc_source_change/(double)(augmented_moves_proposed[3]);
  double swap_prob = (double)nacc_swap/(double)(augmented_moves_proposed[4]);
  double import_swap_prob = (double)nacc_import_swap/(double)(augmented_moves_proposed[5]);
  Rcout << "Beta acceptance = " << beta_prob << ", Lambda acceptance = " << lambda_prob << std::endl;
  Rcout << "Move prob = " << move_prob << ", add prob = " << add_prob << ", remove prob = " << remove_prob 
        << ", source change prob = " << source_change_prob << ", swap prob = " << swap_prob 
        << ", import swap prob = " << import_swap_prob << std::endl;
  */
  
  
}





// [[Rcpp::export]]
void MCMC_Epi_NS(List MCMC_options, IntegerVector t_a, IntegerVector t_c, IntegerVector t_d, IntegerVector source, IntegerMatrix screening_matrix)
{
  // Initialise data
  List data = InitialiseData(t_a, t_c, t_d, source, screening_matrix);
  
  // Load MCMC options
  const int iterations = MCMC_options["iterations"];
  const int population_size = t_a.length();
  const int num_updates = MCMC_options["num_updates"];
  NumericVector prior_parameters = MCMC_options["prior_parameters"];
  NumericVector initial_chain_state = MCMC_options["initial_chain_state"];
  NumericVector parameters = clone(initial_chain_state);
  NumericVector proposal_variance = MCMC_options["proposal_variance"];
  IntegerVector debug_flags = MCMC_options["debug_flags"];
  std::string output_file = MCMC_options["output_file"];
  std::string coltime_file = MCMC_options["coltime_file"];
  std::string source_file = MCMC_options["source_file"];
  std::string loglik_file = MCMC_options["loglik_file"];
  
  // Calculate other variables
  int true_positives = CalculateTruePositives(data);
  int nacc_beta = 0;
  int nacc_add = 0;
  int nacc_move = 0;
  int nacc_remove = 0;
  double loglik = CalculateLogLikelihood_NS(data, parameters);
  IntegerVector Va;
  IntegerVector Vs = CalculateVs(data);
  IntegerVector Vq = CalculateVq(data);
  
  // Write to files
  // Output file
  remove(output_file.c_str());
  std::ofstream myfile; // Define output stream
  myfile.open(output_file.c_str()); // Open file
  assert(myfile.is_open());
  PrintNumVectorToFile(myfile, parameters);
  
  // Colonisation time sum file
  remove(coltime_file.c_str());
  std::ofstream myfile2;
  myfile2.open(coltime_file.c_str());
  assert(myfile2.is_open());
  PrintColonisationTimeSumToFile(myfile2, data);
  
  // Source file
  remove(source_file.c_str());
  std::ofstream myfile3;
  myfile3.open(source_file.c_str());
  assert(myfile3.is_open());
  PrintIntVectorToFile(myfile3, data["source"]);
  
  // loglik file
  remove(loglik_file.c_str());
  std::ofstream myfile4;
  myfile4.open(loglik_file.c_str());
  assert(myfile4.is_open());
  myfile4 << loglik << std::endl;
  
  
  if(!CheckData(data)) stop("inconsistency in data before MCMC");
  
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
      int import_sum = CalculateImportSum(data);
      parameters[1] = R::rbeta(prior_parameters[2] + import_sum, prior_parameters[3] + population_size - import_sum);
    }
    
    loglik = CalculateLogLikelihood_NS(data, parameters);
    
    // Update beta by metropolis hastings
    if(debug_flags[2]==0)
    {
      UpdateTransmissionRate_NS(data, parameters, loglik, proposal_variance[0], prior_parameters[4], nacc_beta);
    }
    
    // Augmented data updates
    if(debug_flags[4]==0)
    {
      for(int update_counter = 0; update_counter < num_updates; update_counter++)
      {
        double w = 0.3;
        double move = R::runif(0.0,3.0);
        if(floor(move) < 1)
        {
          // Move a colonisation time
          MoveColonisationTime_NS(data, parameters, loglik, w, Vq, Va, nacc_move);
        }
        else if(floor(move) < 2)
        {
          // Add a colonsation time
          AddColonisationTime_NS(data, parameters, loglik, w, Vs, Va, nacc_add);
        }
        else if(floor(move) < 3)
        {
          // Remove a colonisation time
          RemoveColonisationTime_NS(data, parameters, loglik, w, Vs, Va, nacc_remove);
        }
        else if(floor(move) < 4)
        {
          // change a source without changing a colonisation time
          ChangeSource_NS(data, parameters, loglik, Vq, Va);
        }
        
        // Write colonisation time sum to file
        PrintColonisationTimeSumToFile(myfile2, data);
        PrintIntVectorToFile(myfile3, data["source"]);
        
        myfile4 << loglik << std::endl;
      }
    }
    
    
    
    // Write parameters to file
    PrintNumVectorToFile(myfile, parameters);
    
    if(i%100==0)
    {
      Rcout << i << std::endl;
    }
    
  }
  
  // Close files
  myfile.close();
  myfile2.close();
  myfile3.close();
  myfile4.close();
  
  double prob = (double)nacc_beta/(double)iterations;
  double move_prob = (double)nacc_move/(double)(iterations * num_updates);
  double add_prob = (double)nacc_add/(double)(iterations * num_updates);
  double remove_prob = (double)nacc_remove/(double)(iterations * num_updates);
  Rcout << "Beta acceptance = " << prob << std::endl;
  Rcout << "Move prob = " << move_prob << ", add prob = " << add_prob << ", remove prob = " << remove_prob << std::endl;
  
  
  
}






