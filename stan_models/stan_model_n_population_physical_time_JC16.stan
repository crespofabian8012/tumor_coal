functions{
  vector rep_vector_times(vector x, int[] times) {
    int N = rows(x);
    int Ntimes = sum(times);
    int sizetimes = size(times);

    vector[Ntimes] y;
    int pos = 1;
    for (n in 1:sizetimes) {
   
      y[pos:(pos+times[n]-1)]= rep_vector(x[n],times[n]);
      pos= pos+times[n];
   
    }
    return y;
  }
vector rep_array_times(real[] x, int[] times) {
    int N = size(x);
    int Ntimes = sum(times);
    int sizetimes = size(times);

    vector[Ntimes] y;
    int pos = 1;
    for (n in 1:sizetimes) {
  
      y[pos:(pos+times[n]-1)]= rep_vector(x[n],times[n]);
      pos= pos+times[n];
   
    }
    return y;
  }
  real oneOnX_log(real x){
		return -log(x);
	}
  int count_valid_elements_vector(vector values){
    int count =1;
    int result;
    if (rows(values)==0){
       result=0;
    }
    else{
       while(count>=1 && count<= rows(values) && !is_nan(values[count])){
        count = count +1;
        }
     result = count -1;
     }
    return(result);
  }
  int find_integer(int x, int[] values){
     int i=1;
     int result= -1;
     while (i <= num_elements(values) && (abs(values[i] -x) >0.00001))
     {
          i=i+1;
     }
     if ((num_elements(values)==0) ||  (i>num_elements(values)))
        result=-1;
      else
         result=i;
      return(result);
  }
  int find(real x, vector values){
     int i=1;
     int result= -1;
     while (i <= rows(values) && (fabs(values[i] -x) >0.00001))
     {
          i=i+1;
     }
     if ((rows(values)==0) ||  (i>rows(values)))
        result=-1;
      else
         result=i;
      return(result);
  }
  int find_sorted(real x, vector sorted_values){
     int i=1;
     int result= -1;
     while (i <= rows(sorted_values)  && (fabs(sorted_values[i] -x) >0.00001) && sorted_values[i]< x)
     {
          i=i+1;
     }
     if ((rows(sorted_values)==0) ||  (i>rows(sorted_values)) || (i>= 1 && i<=rows(sorted_values) && sorted_values[i]> x) )
        result=-1;
      else
         result=i;
      return(result);
  }
  int[]  cumulative_sum_integer(int[] values, int N){
       int result[N];
       int current_sum=0;
       for(i in 1:N) {
            current_sum+=values[i];
            result[i]=current_sum;
       }
       return(result);
  }
  vector  difference_sorted_lists(vector sorted_list, vector sorted_sublist){
       vector[0] result;
       int idx_list=1;
       int idx_sublist=1;
       while(idx_list <= rows(sorted_list)) {
            if (sorted_list[idx_list]< sorted_sublist[idx_sublist]){
              result =append_row(result,sorted_list[idx_list] );
              idx_list= idx_list+1;
            }
            else if (sorted_list[idx_list]> sorted_sublist[idx_sublist]){
              result =append_row(result,sorted_list[idx_list] );
              idx_sublist= idx_sublist+1;
            }
            else{
              idx_list= idx_list+1;
              idx_sublist= idx_sublist+1;
            }
       }
       return(result);
  }
  real model_time_to_standard_time(real t,
         real torigin,
         real delta,
         real K)
         {
          real model_time=0.0;
          real a = exp(delta * t) - 1.0;
          real b = 1.0 - exp(-1.0 * delta * torigin);
          real c = 1.0 - exp(-1.0 * delta * (torigin - t));
          real d;
          real numerator;
          real  denominator;
          if (K==0){
            model_time = a * b / (delta * c);
            return(2.0*model_time/ delta);
          }
         else{
           c = exp(-1.0 * delta * torigin);
           d = 1.0 -c;
           a = ((K / delta) * d) - c;
           b = 1.0 - (K/ delta)*d;
           numerator =(a*exp(delta*t)+b);
           denominator = 1-exp(-delta*(torigin-t));
           numerator =(a*exp(delta*t)+b);
           denominator = 1-exp(-delta*(torigin-t));

           model_time = (2.0 / K )*log(numerator/denominator);
           return(model_time);
            }
         }
    real standard_time_to_model_time(real V,
         real torigin,
         real delta,
         real K)
         {
        real a;
        real b;
        real c;
        real d;
        real e;
        real x;
        real firstTerm;
        real secondTerm = 1.0 / delta;
        real thirdTerm;
        real numerator;
        real denominator;
        real StandardTimeG;
        if (K==0){
           a =  exp(-1.0 * delta * torigin);
           b = (1 - a) * (1 - a) * (1.0 / a);
           c = 1 - a;
           d = (1 - a) * (1.0 / a);
           e = V + d;
           thirdTerm = log(1 - b / e);
           thirdTerm = log(1 - ((1 - a) * (1 - a) * (1.0 / a)) / (V * delta + (1 - a) * (1.0 / a)));
           thirdTerm = log(1 + delta * V - a) - log(1 + (delta * V - 1) * a);
            StandardTimeG = secondTerm * thirdTerm;

        if ( (1 + delta * V - a) <= 0 ||   (1 + (delta * V - 1)*a ) <= 0 ) // do approximation if required
        {
            StandardTimeG = 0.0;
            firstTerm = 0.0;
            secondTerm = 0.0;
            thirdTerm = 0.0;
            a = 0.0;
            b = 0.0;
            c = 0.0;
            d = 0.0;
            e = 0.0;
            a = 1 / delta;
            b = log(1 + delta * V);
            firstTerm = a * b;
            d = (V * V * delta * exp(-1.0 * delta * torigin)) / (1 + V);
            secondTerm =  d;
            StandardTimeG = firstTerm - secondTerm;
           }
       }
       else {
          x = exp(K*V /2.0);
          c = exp(-1.0 * delta * torigin);
          d= 1.0 -c;
          a = ((K/delta) * d ) - c;
          b = 1.0 - (K/ delta)*d;
          numerator = x -b;
          denominator = (x*c) +a;
          thirdTerm= log(numerator/denominator);
          StandardTimeG = secondTerm * thirdTerm;
          }
        return StandardTimeG;
         }
      vector  from_model_time_to_Kingman_coalescent_time(vector coal_times_model_time, real torigin,
         real delta,
         real K){

        vector [rows(coal_times_model_time)] result;
        for(i in 1:rows(coal_times_model_time)){
           if ( !is_nan(coal_times_model_time[i]))
               result[i]= model_time_to_standard_time(coal_times_model_time[i], torigin, delta, K );
           else
              result[i]=-1;
        }
        return(result);
      }
  real log_lik_no_coalescent_between_times(real from,  real to, int number_alive_ind, real torigin, real delta, real  K){
    real  result =  0.0;
    real cumulative_to;
    real cumulative_from;
    if (number_alive_ind>1){
      cumulative_to = model_time_to_standard_time(to, torigin, delta, K);
      cumulative_from = model_time_to_standard_time(from, torigin, delta, K);
      result= log(choose(number_alive_ind,2)) -1.0 *choose(number_alive_ind,2)*(cumulative_to-cumulative_from);
    }
    return(result);
  }
  real conditionalDensityTOrigin_rng(real delta, int sample_size){
        real number_ancestors_population_when_sample_minus_one;
        real U;
        real torigin;
        real torigin_in_model_time_oldest_population;
        real term;
        number_ancestors_population_when_sample_minus_one = 2*sample_size;
        U = gamma_rng(number_ancestors_population_when_sample_minus_one+1, 1);
        term = 1.0- delta /(U+delta);
        if (term <= 0){
          reject("term must not be negative; found term=", term);
        }
        else{
            torigin= -(1.0 / delta)*log(term);
        }
        return(torigin);
  }
  real conditionalDensityTOrigin_pdf( real y, real delta,  int sample_size) {
      real result;
      real term1;
      real term2;
      real term3;
      real partial;
      term1 = exp(-1.0*delta*y);
      term2 = delta * term1;
      term3 = 1.0-term1;
      partial = delta * term2 /(term3 * term3);
      partial = partial * exp( -1 * term2/term3);

      result  = partial;
      return(result);
   }
  real conditionalDensityTOrigin_lpdf( real y, real delta,  int sample_size) {
      real result;
      result = log(conditionalDensityTOrigin_pdf(y, delta,   sample_size));
      return(result);
   }
  real conditionalDensityTOrigin_cdf(real y, real delta,  int sample_size) {
     real term1;
     real term2;
     real term3;
     real result;
     term1 = exp(-1.0 * delta * y);
     term2 = delta * term1;
     term3 = 1.0-term1;
     result = exp( -1.0 * term2/term3);

     return(result);
  }
    real conditionalDensityTOrigin_lcdf(real y, real delta,  int sample_size) {
     real term1;
     real term2;
     real term3;
     real result;
     term1 = exp(-1.0 * delta * y);
     term2 = delta * term1;
     term3 = 1.0-term1;
     result =  -1.0* term2/term3;
    return(log(result));
 }
  real conditionalDensityTOrigin_lccdf(real y, real delta,  int sample_size) {
    real term1;
    real term2;
    real term3;
    real result;
    result = conditionalDensityTOrigin_cdf(y, delta, sample_size );
    return(log(1-result));
 }
  real log_h(real t,
         real torigin,
         real delta,
         real K){

    real a = 1.0 - exp(-1.0 * delta * (torigin - t));
    real first_term = 2.0 * log(a);
    real second_term = -1.0 * delta * t;
    real third_term = exp(delta * t);
    real above_term = first_term + second_term;
    real b = 1.0 - exp(-1.0 * delta * torigin);
    real  below_term = 2.0 * log(b);
    real extra_term = log(1+ (K /delta)*b*(third_term -1.0) / a );
    real logH;
    if (K==0)
        logH = above_term - below_term;
    else
        logH =  above_term - below_term + extra_term;
    return(logH);
  }
  real log_prob_father_population(int order, int index_oldest_population, real father_torigin, real father_delta,
                                 real x_father, real torigin, real delta, real x, vector deltas,
                                 vector torigins_in_model_time,
                                vector torigins_in_model_time_oldest_population_par, real K, int N, vector pop_sizes_proportion){
                                 // , simplex pop_sizes_proportion
         real log_numerator=0.0;
         real log_denominator=0.0;
         real logh;
         real temp;
         real result;
    
         if (order != N ){

                for(i in (order+1):N){
                  int pos = i;
                  log_denominator = log_denominator + log( pop_sizes_proportion[pos]);
                  temp= torigin * x / pop_sizes_proportion[pos];
                  logh=log_h(temp, torigins_in_model_time[pos], deltas[pos], K);
                  log_denominator = log_denominator  + temp;
   
                }
         }

         log_numerator = log_numerator + log( x_father);
         temp= torigin * x / x_father;
         logh=log_h(temp, father_torigin, father_delta, K);
         log_numerator = log_numerator  + temp;
         

         result= log_numerator - log_denominator;

        return(result);
  }
   real log_likelihood_subpopulation(int order, int sample_size, real K, real delta, real torigin, int[] positions,
                                     vector sorted_coal_time_with_previous_torigins_in_model_time, vector child_torigins_in_model_time){

    real log_lik = 0.0;
    real current_time;
    int  alive_cells;
    real term_only_after_first_coal_event;
    real zero =0.0;
    int currentCoalescentEventInThisEpoch=-1;
    vector[rows(sorted_coal_time_with_previous_torigins_in_model_time)+1] padded_coalescent_times =
                                                  append_row(zero, sorted_coal_time_with_previous_torigins_in_model_time);
    real last_event_time_before_migration=0.0;


    alive_cells = sample_size;//alive_cells goes from sample_size downto 2
    for(j in 1:rows(padded_coalescent_times)){
      current_time = padded_coalescent_times[j];

      //print("current_time=",current_time);

      if (fabs(current_time - torigin)<0.00000001)
        return log_lik;

      if (current_time < 0.0) continue;//those are not real inmigrants

      if (rows(child_torigins_in_model_time) > 0 && find(current_time, child_torigins_in_model_time)!= -1){//= -1
        

        if (alive_cells >1){
            log_lik+= log_lik_no_coalescent_between_times( last_event_time_before_migration,  current_time,  alive_cells,  torigin,  delta,  K);
            alive_cells= alive_cells+1;
            currentCoalescentEventInThisEpoch=1;
            last_event_time_before_migration= current_time;
            }
      }
      else{

         //alive_cells =(sample_size-j+1);//alive_cells goes from sample_size downto 2
        if (alive_cells >1) {
            log_lik = log_lik + log(alive_cells*(alive_cells-1)/2.0);
            log_lik = log_lik + log(2.0)-log_h(current_time,torigin, delta, K);
            if (fabs(current_time)<0.00001){
              term_only_after_first_coal_event=0.0;
            }
            else{
                   term_only_after_first_coal_event= model_time_to_standard_time(last_event_time_before_migration,torigin, delta, K);
              }
             log_lik = log_lik - (alive_cells*(alive_cells-1)/2.0)*(model_time_to_standard_time(current_time,torigin, delta, K)-
                  term_only_after_first_coal_event);
              alive_cells= alive_cells-1;
              currentCoalescentEventInThisEpoch= currentCoalescentEventInThisEpoch+1;
              last_event_time_before_migration= current_time;
            }
         }
      }

      return(log_lik);
    }
    vector structured_coalescent_rng(vector deltas,vector torigins_in_model_time, vector torigins_in_model_time_oldest_population,
                                   vector pop_sizes_proportion, int[] sample_sizes, int[] number_internal_nodes,
                                   int[] indexes_father_populations, int[] pos_parent_child_pop_MRCA_in_parent_pop, int[] cum_internal_nodes, 
                                   int[] number_child_pops,
                                   int[] child_pop_indexes,
                                   real K ){
    int total_sample_size=sum(sample_sizes);
    int N = rows(deltas);
    int current_pos=1;
    int current_pos_migrations=1;
    int number_inmigrants=1;
    int index_next_immigrants=1;
    int global_index_next_immigrants=1;
    //counting torigin as immigrant
    int  indexes_inmigrants = 1;
    vector[total_sample_size-1] coal_event_times_model_time_oldest_population;

    coal_event_times_model_time_oldest_population=rep_vector(torigins_in_model_time_oldest_population[1]+1,total_sample_size-1);

    print("inside structured_coalescent_rng=",sample_sizes);
    for(i in 1:num_elements(sample_sizes)){
        int sample_size= sample_sizes[i];
        int current_sample_size= sample_size;
        int index_father_population=0;
        current_pos=1;
        if (i==1){
            vector[sample_size-1] number_ancestors_sample;
            vector[sample_size-1] rate_exp;
            vector[sample_size-1] times;
            //vector[sample_size-1] reverse_times;
            vector[sample_size-1] cum_sum_times;

            for(j in 2:sample_size){
               number_ancestors_sample[j-1]=j;
               rate_exp[j-1]=0.5*j*(j-1);
               times[j-1] = exponential_rng(rate_exp[j-1]);
            }
            //print("times=",times);
            for(j in 1:(sample_size-1)){
              cum_sum_times[j] = sum(times[j:(sample_size-1)]);
              coal_event_times_model_time_oldest_population[ sample_size-j]= standard_time_to_model_time(cum_sum_times[j],
                                                                                                 torigins_in_model_time[i],
                                                                                                 deltas[i],
                                                                                                 K);

            }
            //print("i=1 coal_event_times_model_time=", coal_event_times_model_time_oldest_population);
            coal_event_times_model_time_oldest_population[1:(sample_size-1)] = (pop_sizes_proportion[i]/pop_sizes_proportion[N]) *
                                                     (coal_event_times_model_time_oldest_population[1:(sample_size-1)]);
            coal_event_times_model_time_oldest_population[sample_size]= torigins_in_model_time_oldest_population[i];
            current_pos = sample_size+1;
            //print("coal_event_times_model_time_oldest_population=", coal_event_times_model_time_oldest_population);
            index_father_population=2;
        }
    }
    return(coal_event_times_model_time_oldest_population);
  }
  real structured_coalescent_lpdf(vector coal_event_times_model_time, vector deltas, vector torigins_in_model_time,
                              vector torigins_in_physical_time,  vector pop_sizes_proportion, int[] sample_sizes, int[] number_internal_nodes,
                               int[] indexes_father_populations, int[] pos_parent_child_pop_MRCA_in_parent_pop, int[] cum_internal_nodes, 
                               int[] number_child_pops,
                               int[] child_pop_indexes,
                               real K){

    real log_lik=0.0;
    int cum_sample_size=0;
    int cum_number_coal=0;
    int first_position_coal_event=1;
    int N= rows(deltas);
    int pos_next_torigin = 0;
    int cum_sum_number_child_pops[N] = cumulative_sum_integer(number_child_pops, N);
    int k=1;
    int pos_idx_child_pop=1;

    while(cum_sum_number_child_pops[k]==0){
      k=k+1;

    }


    cum_sample_size=0;

    for(i in 1:num_elements(sample_sizes)){
         int order=i;
         int index_father_population = indexes_father_populations[i];
         int number_events= number_internal_nodes[i];
         int positions[sample_sizes[i]] = rep_array(1,sample_sizes[i] );
         vector[number_events+number_child_pops[i]] current_coal_event_times_with_immigrants;
         vector[number_events+number_child_pops[i]] sorted_current_coal_event_times_with_immigrants;
         vector[number_child_pops[i]] child_pops_inmigration ;
         vector[number_child_pops[i]] sorted_child_pops_inmigration ;
         int start_idx=1;
         int start_idx_number_pop=1;
         pos_idx_child_pop=1;
         
         //print("k=",k);
         if (i>=2){
          start_idx = cum_internal_nodes[i-1]+1;
          start_idx_number_pop= k;
         }
  

         current_coal_event_times_with_immigrants[1:number_events] = (coal_event_times_model_time[ start_idx:(cum_internal_nodes[i])]) ;

         //print("current_coal_event_times_with_immigrants[1:number_events] ", current_coal_event_times_with_immigrants[1:number_events]);

         if (number_child_pops[i]>0){

            
            for(j in 1:num_elements(indexes_father_populations)){

               if(indexes_father_populations[j]==i){


                  
                   child_pops_inmigration[pos_idx_child_pop]= (1.0/pop_sizes_proportion[i])*torigins_in_physical_time[j];
                   
                    //print(" child_pops_inmigration[pos_idx_child_pop] con physical time",child_pops_inmigration[pos_idx_child_pop]);
                   child_pops_inmigration[pos_idx_child_pop]= (pop_sizes_proportion[j]/pop_sizes_proportion[i])*torigins_in_model_time[j];

       
                     
                     if(child_pops_inmigration[pos_idx_child_pop]> current_coal_event_times_with_immigrants[pos_parent_child_pop_MRCA_in_parent_pop[j]]){

                    //print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!wrong position of time of origin of child pop");
                      }
                    //}
                 // }
                  

                   pos_idx_child_pop=pos_idx_child_pop+1;
               }

            }
           
     
         }

           current_coal_event_times_with_immigrants[(number_events+1):(number_events+number_child_pops[i])]=child_pops_inmigration;
           
           sorted_current_coal_event_times_with_immigrants= sort_asc(current_coal_event_times_with_immigrants);
           sorted_child_pops_inmigration  = sort_asc(child_pops_inmigration);
       

           log_lik= log_lik+ log_likelihood_subpopulation(i, sample_sizes[i], K, deltas[i], torigins_in_model_time[i], positions,
                                      sorted_current_coal_event_times_with_immigrants,
                                      sorted_child_pops_inmigration
                                      //rep_vector(1,1)
                                      );

             if (i<(N-1) && index_father_population>0 && N>2){
              log_lik+= log_prob_father_population(order, N, torigins_in_model_time[indexes_father_populations[i]], deltas[indexes_father_populations[i]],
                                           pop_sizes_proportion[indexes_father_populations[i]],
                                           torigins_in_model_time[i], deltas[i], pop_sizes_proportion[i], deltas,  torigins_in_model_time, torigins_in_physical_time,K,N,  pop_sizes_proportion);
            }

    }
    return(log_lik);
  }
}
data{
  real<lower=0.0, upper=2.0> K;//parameter for the model M_K
  int<lower=2> N;//number of populations
  int<lower=2> total_sample_size;
  int<lower=0> L;// alignment length
  int<lower=0> ids_MRCA[N];
  int<lower=0> sample_sizes[N];
  int<lower=1> number_internal_nodes[N];
  int indexes_father_populations[N];
  //for example [0,1,1,2] means the population 2 and 3 have father population 1(oldest population from present time)
  //and population 4 has parent population 2
  int number_child_pops[N];
  int child_pop_indexes[sum(number_child_pops)];
  int  pos_parent_child_pop_MRCA_in_parent_pop[N];

  int<lower=0,upper=16> genotype_tipdata[total_sample_size,L];
  int<lower=1> site_weights[L];
  int<lower=0,upper=2*total_sample_size> topology[total_sample_size-1,3];
  int<lower=0,upper=total_sample_size> tip_association[total_sample_size];
  matrix[2*total_sample_size-2, total_sample_size-1 ] coal_times_to_branch_lengths;
}
transformed data{
    int number_branches = 2*total_sample_size-2;
    int number_of_tips_below[2*total_sample_size-1];
    
    int number_of_coalescent_times_below[2*total_sample_size-1] ;
    vector[16] genotype_tip_partials[2*total_sample_size,L];
    int tip_index;
    

    matrix[total_sample_size-1, 2*total_sample_size-1] indexes_nodes_below= rep_matrix(0, total_sample_size-1,2*total_sample_size-1 );
    //int map_internal_node_topology_row[2*total_sample_size-1]= rep_array(0,2*total_sample_size-1);
    //matrix[number_branches, total_sample_size-1 ] coal_times_to_branch_lengths= rep_matrix(0, number_branches, total_sample_size-1);
    vector [3*total_sample_size-4] w;
    int v[3*total_sample_size-4];
    int u[number_branches+1];
    int max_sample_size= max(sample_sizes);
    int cum_sample_sizes[N];
    int cum_internal_nodes[N];
    int cum_number_child_pops[N];

    cum_sample_sizes= cumulative_sum_integer(sample_sizes, N);
    cum_internal_nodes= cumulative_sum_integer(number_internal_nodes, N);
cum_number_child_pops= cumulative_sum_integer(number_child_pops, N);
 
    for(i in 1:(total_sample_size)){
       number_of_tips_below[i]=1;
       number_of_coalescent_times_below[i]=0;
    }
      for(i in 1:(total_sample_size-1)){
        int parent_node_idx = topology[i,1];
        int left_node_idx = topology[i,2];
        int right_node_idx = topology[i,3];
        //int current_row = i;
        number_of_tips_below[parent_node_idx] = number_of_tips_below[left_node_idx]+ number_of_tips_below[right_node_idx];
        number_of_coalescent_times_below[parent_node_idx] = number_of_coalescent_times_below[left_node_idx]+ number_of_coalescent_times_below[right_node_idx]+1;
     
}

    w=csr_extract_w(coal_times_to_branch_lengths);
    v=csr_extract_v(coal_times_to_branch_lengths);
    u=csr_extract_u(coal_times_to_branch_lengths);

    for( n in 1:total_sample_size ) {//rows
          for( i in 1:L ) {//columns
                 genotype_tip_partials[n,i]= rep_vector(0.0, 16);
                 tip_index = tip_association[n];
                 if (genotype_tipdata[tip_index,i]!= 0){
                   genotype_tip_partials[n,i][genotype_tipdata[tip_index,i] ] = 1.0;
                 }
                 else{
                   genotype_tip_partials[n,i] = rep_vector(1.0, 16);
                 } 
          }
      }

    

}
parameters{
  real<lower=0.0> hyper_parameters_exponential[N];
  //real<lower=0.0> hyper_parameters_exponential_theta;
  vector<lower=0.001>[N] gammas;
  //the first gamma corrsespond to the oldest population
  
  real<lower=0.0> theta;
  real<lower=0.0>  torigin_oldest_population;
  
  simplex[N] pop_sizes_proportion;
  simplex[total_sample_size] population_simplexes[N];
  vector<lower=0, upper=1>[N-1] positions_torigins_in_intervals; 
  
  //vector<lower=0.0>[total_sample_size-1+N] gamma_variables;
}
transformed parameters {
  //declarations
 
  vector[total_sample_size-1] coal_event_times_physical_time;
  vector[total_sample_size-1] coal_event_times_model_time;
  vector<lower=0.0>[N] torigins_in_model_time;
  vector<lower=0.0>[N] torigins_in_physical_time;
  vector<lower=0>[sum(number_internal_nodes)+N] pi;

  torigins_in_model_time[N] = torigin_oldest_population;
  torigins_in_physical_time[N] = torigins_in_model_time[N]* pop_sizes_proportion[N];
  
 
   for(i in 1:N){
       int start_idx=1;
       int end_idx=1;
       int no_internals ;
        if (i>=2){
          start_idx = cum_sample_sizes[i-1]+1+cum_number_child_pops[i-1];
          
        }
        end_idx = cum_sample_sizes[i]+cum_number_child_pops[i];
    
      
       no_internals = number_internal_nodes[i];


       pi[start_idx:end_idx] = append_row(head(population_simplexes[i], no_internals), sum(tail(population_simplexes[i], total_sample_size - no_internals  )));
    
   }


   for (k in 1:(N-1)) {
      int j= N-k; 
      int idx_parent_pop = indexes_father_populations[j];
      vector [number_internal_nodes[idx_parent_pop]+1] pi_father_pop;
      vector[pos_parent_child_pop_MRCA_in_parent_pop[j]] head_cum_sum;
      real lower_coal_time=0;
      real upper_coal_time;
      real unif;
      int end_idx=1;
      int no_internals = number_internal_nodes[idx_parent_pop-1];
      //real torigins_in_model_time_parent_population;
      int start_idx=1;
        if (idx_parent_pop>=2){
          start_idx = cum_sample_sizes[idx_parent_pop-1]+1+cum_number_child_pops[idx_parent_pop-1];
        }
       end_idx = cum_sample_sizes[idx_parent_pop]+cum_number_child_pops[idx_parent_pop];

      pi_father_pop = pi[start_idx:end_idx];
      
   
      head_cum_sum= head(cumulative_sum(pi_father_pop), pos_parent_child_pop_MRCA_in_parent_pop[j]);
 
      
      if (pos_parent_child_pop_MRCA_in_parent_pop[j]==1){
         lower_coal_time = 0;

      }else{
        lower_coal_time = 0;// head_cum_sum[pos_parent_child_pop_MRCA_in_parent_pop[j]-1]; 
      }
   

      upper_coal_time = head_cum_sum[pos_parent_child_pop_MRCA_in_parent_pop[j]]; 
      unif= positions_torigins_in_intervals[j];

      torigins_in_physical_time[j] =  torigins_in_physical_time[idx_parent_pop] *unif*upper_coal_time ;
  
    
      torigins_in_model_time[j] = torigins_in_physical_time[j] / pop_sizes_proportion[j];

   }
 

   for (j in 1:N) {
        int start_idx=1;
        int end_idx=1;
        int start_coal_idx=1;
        int end_coal_idx=1;
        if (j>=2){
          start_idx =cum_sample_sizes[j-1]+1+cum_number_child_pops[j-1];
          start_coal_idx= cum_internal_nodes[j-1]+1;
        }

    
        end_idx = cum_sample_sizes[j]+cum_number_child_pops[j];

        end_coal_idx= cum_internal_nodes[j];

       
       coal_event_times_model_time[start_coal_idx:end_coal_idx] = head(cumulative_sum(pi[start_idx:end_idx]), number_internal_nodes[j]) *  torigins_in_model_time[j];
       
       coal_event_times_physical_time[start_coal_idx:end_coal_idx] = head(cumulative_sum(pi[start_idx:end_idx]), number_internal_nodes[j]) *  torigins_in_physical_time[j]  ;
     
    }

}
model{
  real a;
  real a_divide_by_3;
  vector[16] left;
  vector[16] right;
  vector[16] partials[total_sample_size,L];  // partial probabilities for the S tips and S-1 internal nodes
  matrix[16,16] p_matrices[number_branches]; // finite-time transition matrices for each branch
  vector[number_branches] branch_lengths;
  vector[total_sample_size-1] J_diag;
  vector[number_branches] J_diag_theta;
  vector[N] tot_pop_size= rep_vector(1000.0, N);
  int pos = 1;
    
  hyper_parameters_exponential~gamma(rep_vector(0.001,N), rep_vector(0.001,N));
  gammas ~ exponential(hyper_parameters_exponential);
  
  theta ~ exponential(1);

   pop_sizes_proportion~ dirichlet(tot_pop_size .* to_vector(sample_sizes));

  
  torigin_oldest_population~conditionalDensityTOrigin(gammas[N], sample_sizes[N]);

   
   for (j in 1:N) {
     population_simplexes[j]~ dirichlet(append_row(rep_vector(1, number_internal_nodes[j] ), 
                     rep_vector(1.0/(total_sample_size - number_internal_nodes[j]), total_sample_size  - number_internal_nodes[j])));
     
   }

  J_diag = rep_vector_times(torigins_in_model_time, number_internal_nodes);
  //this is a Jacobian from the simplexes to the  coalescent event times 

  target += sum(log(J_diag));

 for (j in 1:(N-1)) {
    int index_father_pop = indexes_father_populations[j];
    int number_internal_nodes_parent_pop = number_internal_nodes[index_father_pop];
    vector[pos_parent_child_pop_MRCA_in_parent_pop[j]] head_cum_sum;
    //vector [number_internal_nodes[index_father_pop]] pi_father_pop = pi[start_idx:cum_internal_nodes[index_father_pop]];
    
    vector[number_internal_nodes_parent_pop] father_pop_simplex = population_simplexes[index_father_pop][1:number_internal_nodes_parent_pop];
    head_cum_sum= head(cumulative_sum(father_pop_simplex), pos_parent_child_pop_MRCA_in_parent_pop[j]);

    //print("father_pop_simplex=",father_pop_simplex);
    //this if ro the change from torigins_in_model_time[N] to torigins_in_model_time_parent_population[N]
    //target +=  log(pop_sizes_proportion[N]);

    target += log(torigins_in_model_time[index_father_pop])+log( head_cum_sum[pos_parent_child_pop_MRCA_in_parent_pop[j]]);
    //+log(torigins_in_model_time[index_father_pop]);
    //jacobian of the transformation from u to times origin child populations in the same time unit of the father T
    target +=  log(pop_sizes_proportion[index_father_pop])-log(pop_sizes_proportion[j]);
    //jacobian of the transformation to convert to the unit of time origin in the child population

    torigins_in_model_time[j]~conditionalDensityTOrigin(gammas[j], sample_sizes[j]);
}



coal_event_times_model_time~structured_coalescent(gammas, torigins_in_model_time,
                    torigins_in_physical_time,
                    pop_sizes_proportion, sample_sizes, number_internal_nodes, indexes_father_populations, pos_parent_child_pop_MRCA_in_parent_pop, cum_internal_nodes, 
                    number_child_pops,
                    child_pop_indexes,
                    K);

branch_lengths =   csr_matrix_times_vector(number_branches, total_sample_size-1, w,  v, u, coal_event_times_physical_time);

 for( b in 1:number_branches){
      a = (1.0/16.0)*(1-exp(-16.0*branch_lengths[b]*theta/15.0));
      p_matrices[b]= rep_matrix(a, 16, 16);
      for(i in 1:16){
        p_matrices[b][i,i]= 1.0-15*a;
      }
    }

for( i in 1:L ){//columns
   
		for( n in 1:(total_sample_size-1) ) {//rows

          left  = (topology[n,2] > total_sample_size)? partials[topology[n,2]-total_sample_size,i] :  genotype_tip_partials[topology[n,2],i];
          right = (topology[n,3] > total_sample_size)? partials[topology[n,3]-total_sample_size,i] :  genotype_tip_partials[topology[n,3],i];
          
          partials[topology[n,1]-total_sample_size,i] = (p_matrices[2*n-1]*left) .* (p_matrices[2*n]*right);
        }
        for(j in 1:16){

        	partials[2*total_sample_size-total_sample_size,i][j] = partials[topology[total_sample_size-1,1]-total_sample_size,i][j] * (1.0/16);
        }
        target += site_weights[i]*log(sum(partials[total_sample_size,i]));
   }
}

generated quantities {
  vector<lower=0.001>[N] unscaled_gammas;
 

  for (i in 1:N) {
    unscaled_gammas[i] = gammas[i]/pop_sizes_proportion[i] ;
  }

}
