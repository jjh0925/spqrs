//=========================================================================================
#include <Rcpp.h>
using namespace Rcpp;
//=========================================================================================
// functions used in spqrs
NumericMatrix bsplines(NumericVector x, NumericVector t,
                       int degree, int derivative);
NumericVector bspline(NumericVector x, NumericVector t,
                      int degree, int j, int derivative);
double bsp(double x, NumericVector t, int degree, int j);
double dbsp(double x, NumericVector t, int degree, int j, int derivative);
NumericMatrix jump_bsplines(NumericVector t, int degree);
NumericVector dim2knots(NumericVector predictor, int dimension, int degree);
NumericVector knots2t(NumericVector knots, int degree);
NumericVector lambdas_all(NumericVector response,
                               int number_lambdas,
                               double lambda_max,
                               double epsilon_lambda);
double check(double u, double tau);
double Check(NumericVector v, double tau);
double find_solution(NumericVector v,
                     NumericVector w,
                     double tau,
                     double tau_penalty,
                     double lambda,
                     NumericVector a,
                     NumericVector c,
                     IntegerVector pen_type);
bool find_slope_v(NumericVector v,
                  NumericVector w,
                  double tau,
                  double tau_penalty,
                  double lambda,
                  NumericVector a,
                  NumericVector c,
                  int i, IntegerVector pen_type);
bool find_slope_c(NumericVector v,
                  NumericVector w,
                  double tau,
                  double tau_penalty,
                  double lambda,
                  NumericVector a,
                  NumericVector c,
                  int l, IntegerVector pen_type);
double type_forward(double target, double tau_penalty, double lambda, double a, double c);
double type_backward(double target, double tau_penalty, double lambda, double a, double c);
int find_interval(NumericVector y, double z);
int find_interval_weight(NumericVector w, double tau);
// array operations
ListOf<NumericVector> NumericVectors(int size);
ListOf<IntegerVector> IntegerVectors(int size);
IntegerVector support_of_vector(NumericVector v);
IntegerVector bubble_order(NumericVector vec);
NumericVector subsetNumVec(NumericVector x, IntegerVector index);
IntegerVector subsetIntVec(IntegerVector x, IntegerVector index);
//=========================================================================================
// Main function
// Simultanious Penalized Quantile Regression Splines
// [[Rcpp::export]]
List spqrs(NumericVector response,
           NumericVector predictor,
           NumericVector tau,
           int           degree,
           int           dimension = 10,
           double        tau_penalty = 0.5,
           int           number_lambdas = 3,
           double        lambda_max = 1e+2,
           double        epsilon_lambda = 1e-10,
           int           maxiter = 500,
           double        epsilon_iterations = 1e-6,
           bool          non_crossing = false)
{
   Rcout << "==================================================\n";
   Rcout << "Simultanious Penalized Quantile Regression Splines\n";
   Rcout << "Version 0.1 by Jae-Hwan Jhong (July 04, 2018)\n";
   Rcout << "Department of Statistics, Korea University, Korea\n";
   Rcout << "==================================================\n";
   // dimension
   int sample_size = response.size();
   int order = degree + 1;
   int number_taus = tau.size();
   // int tau_index, lambda_index;
   IntegerVector dimension_vec(number_taus);
   for (int i = 0; i < number_taus; i++)
      dimension_vec[i] = dimension;
   // initial knots
   ListOf<NumericVector> knots_list = NumericVectors(number_taus);
   ListOf<NumericVector> t_list = NumericVectors(number_taus);
   ListOf<NumericVector> coefficients_list = NumericVectors(number_taus);
   ListOf<NumericVector> residuals_list = NumericVectors(number_taus);
   for (int i = 0; i < number_taus; i++)
   {
      knots_list[i] = dim2knots(predictor, dimension_vec[i], degree);
      t_list[i] = knots2t(knots_list[i], degree);
      NumericVector initial_coefficients(dimension_vec[i]);
      coefficients_list[i] = initial_coefficients;
      residuals_list[i] = clone(response);
   }
   // all lambdas
   NumericVector lambdas = lambdas_all(response, number_lambdas, lambda_max, epsilon_lambda);
   // local variables
   List results(number_lambdas + 1);
   NumericVector b_j(sample_size);
   NumericVector v_j(sample_size);
   NumericVector w_j(sample_size);
   NumericVector residuals_j(sample_size);
   // Criterion
   NumericVector aic_vector(number_lambdas);
   NumericVector bic_vector(number_lambdas);
   NumericVector dimension_vector(number_lambdas);
   NumericVector R_vec(number_taus);
   NumericVector R_lambda_vec(number_taus);
   double lambda, R, R_lambda, store_R_lambda, sol_candidate;
   int i, iter, j, k, l, q;
   int number_pruning = 10;
   // initial computations
   IntegerVector number_penalty_vec(number_taus);
   for (i = 0; i < number_taus; i++)
      number_penalty_vec[i] = dimension_vec[i] - order;
   Rcout << "tau = " << tau << std::endl;
   NumericMatrix tilde_B_km1 = bsplines(predictor, knots_list[1], degree, 0);
   NumericMatrix tilde_B_k = bsplines(predictor, knots_list[1], degree, 0);

   // Module FIT
   for (int lambda_index = 0; lambda_index < number_lambdas; lambda_index++)
   {
      Rcout << "\r" << "lambda_index = " << lambda_index;
      lambda = lambdas[lambda_index];
      // fit a spqrs corresponding to the lambda
      for (iter = 0; iter < maxiter; iter++)
      {
         // Module TAUS
         for (int tau_index = 0; tau_index < number_taus; tau_index++)
         {
            NumericMatrix basis = bsplines(predictor, t_list[tau_index], degree, 0);
            NumericMatrix jump = jump_bsplines(t_list[tau_index], degree);
            ListOf<IntegerVector> supports = IntegerVectors(dimension_vec[tau_index]);
            for (j = 0; j < dimension_vec[tau_index]; j++)
               supports[j] = support_of_vector(basis(_, j));
            // Module UPDATE
            if (non_crossing)
            {
               if (tau_index > 0)
               {
                  NumericVector sum_knots = union_(knots_list[tau_index - 1], knots_list[tau_index]);
                  sum_knots = unique(sum_knots);
                  sum_knots.sort();
                  tilde_B_km1 = bsplines(sum_knots, t_list[tau_index - 1], degree, 0);
                  tilde_B_k = bsplines(sum_knots, t_list[tau_index], degree, 0);
               }
            }
            for (j = 0; j < dimension_vec[tau_index]; j++)
            {
               b_j = basis(_, j);
               // partial residual y_ij
               for (k = 0; k < supports[j].size(); k++)
                  residuals_list[tau_index][supports[j][k]] += b_j[supports[j][k]] * coefficients_list[tau_index][j];
               // compute v and w
               w_j = b_j[supports[j]];
               residuals_j = residuals_list[tau_index][supports[j]];
               v_j = residuals_j / w_j;
               IntegerVector order_v = bubble_order(v_j);
               w_j = subsetNumVec(w_j, order_v);
               v_j.sort();
               if (number_penalty_vec[tau_index] == 0)
               {
                  // weighted quantile
                  q = find_interval_weight(w_j, tau[tau_index]);
                  sol_candidate = v_j[q];
                  if (tau_index == 0)
                     coefficients_list[tau_index][j] = sol_candidate;
                  if (!non_crossing)
                     coefficients_list[tau_index][j] = sol_candidate;
               }
               else
               {
                  // compute a, c
                  // compute a, d vector for penalty term
                  int nonzero_size = 0;
                  NumericVector row_jump(jump(j, _));
                  for (l = 0; l < row_jump.size(); l++)
                  {
                     if (std::abs(row_jump[l]) > 1e-6)
                        nonzero_size++;
                  }
                  // memory for a, c
                  NumericVector c(nonzero_size);
                  NumericVector a(nonzero_size);
                  IntegerVector nonzero_index(nonzero_size);
                  int ss = 0;
                  for (l = 0; l < row_jump.size(); l++)
                  {
                     if (std::abs(row_jump[l]) > 1e-6)
                     {
                        a[ss] = row_jump[l];
                        nonzero_index[ss] = l;
                        ss++;
                     }
                  }
                  // compute c
                  for (k = 0; k < nonzero_size; k++)
                     c[k] = coefficients_list[tau_index][j] - sum(jump(_, nonzero_index[k]) * coefficients_list[tau_index]) / a[k];
                  // abs a after compute c
                  a = abs(a);
                  // pen_type
                  IntegerVector pen_type(nonzero_size);
                  for (k = 0; k < nonzero_size; k++)
                     if (jump(j, nonzero_index[k]) > 0)
                        pen_type[k] = 1;
                  sol_candidate = find_solution(v_j, w_j, tau[tau_index], tau_penalty, lambda, a, c, pen_type);
                  if (tau_index == 0)
                     coefficients_list[tau_index][j] = sol_candidate;
                  if (!non_crossing)
                     coefficients_list[tau_index][j] = sol_candidate;
               }
               // Non-Crossing constraint for constant and linear
               if (non_crossing)
               {
                  if (tau_index > 0)
                  {
                     int nonzero_size_nc = 0;
                     NumericVector col_tilde_B_k(tilde_B_k(_, j));
                     for (l = 0; l < col_tilde_B_k.size(); l++)
                     {
                        if (col_tilde_B_k[l] > 1e-6)
                           nonzero_size_nc++;
                     }
                     // memory for a_nc, c_nc
                     NumericVector c_nc(nonzero_size_nc);
                     NumericVector a_nc(nonzero_size_nc);
                     IntegerVector nonzero_nc_index(nonzero_size_nc);
                     int ss = 0;
                     for (l = 0; l < col_tilde_B_k.size(); l++)
                     {
                        if (col_tilde_B_k[l] > 1e-6)
                        {
                           a_nc[ss] = col_tilde_B_k[l];
                           nonzero_nc_index[ss] = l;
                           ss++;
                        }
                     }
                     // compute c_nc
                     for (k = 0; k < nonzero_size_nc; k++)
                     {
                        c_nc[k] = coefficients_list[tau_index][j] +
                           (sum(tilde_B_km1(nonzero_nc_index[k], _) * coefficients_list[tau_index - 1]) -
                           sum(tilde_B_k(nonzero_nc_index[k], _) * coefficients_list[tau_index])) / a_nc[k];
                     }
                     double max_c_nc = max(c_nc);
                     if (sol_candidate < max_c_nc)
                        coefficients_list[tau_index][j] = max_c_nc;
                     else
                        coefficients_list[tau_index][j] = sol_candidate;
                  }
               }
               // recompute residual
               for (k = 0; k < supports[j].size(); k++)
                  residuals_list[tau_index][supports[j][k]] -= b_j[supports[j][k]] * coefficients_list[tau_index][j];
            }
            // Module Prune
            if (iter <= number_pruning)
            {
               if (number_penalty_vec[tau_index] > 0)
               {
                  NumericVector penalty(number_penalty_vec[tau_index]);
                  for (k = 0; k < number_penalty_vec[tau_index]; k++)
                     penalty[k] = sum(jump(_, k) * coefficients_list[tau_index]);
                  IntegerVector penalty_check(number_penalty_vec[tau_index]);
                  penalty_check[abs(penalty) < 1e-6] = 1;
                  if (sum(penalty_check) > 0)
                  {
                     // re-compute
                     IntegerVector prune_coef(dimension_vec[tau_index]);
                     IntegerVector prune_knots(number_penalty_vec[tau_index] + 2);
                     for (k = 0; k < number_penalty_vec[tau_index]; k++)
                     {
                        if (std::abs(penalty[k]) < 1e-6)
                        {
                           prune_coef[k] = 1;
                           prune_knots[k + 1] = 1;
                        }
                     }
                     coefficients_list[tau_index] = coefficients_list[tau_index][prune_coef == 0];
                     knots_list[tau_index] = knots_list[tau_index][prune_knots == 0];
                     dimension_vec[tau_index] = coefficients_list[tau_index].size();
                     // fitted_values
                     NumericVector fitted_values(sample_size);
                     number_penalty_vec[tau_index] = dimension_vec[tau_index] - order;
                     t_list[tau_index] = knots2t(knots_list[tau_index], degree);
                     basis = bsplines(predictor, t_list[tau_index], degree, 0);
                     if (number_penalty_vec[tau_index] > 0)
                        jump = jump_bsplines(t_list[tau_index], degree);
                     for (j = 0; j < dimension_vec[tau_index]; j++)
                     {
                        supports[j] = support_of_vector(basis(_, j));
                        b_j = basis(_, j);
                        for (k = 0; k < supports[j].size(); k++)
                           fitted_values[supports[j][k]] += b_j[supports[j][k]] * coefficients_list[tau_index][j];
                     }
                     residuals_list[tau_index] = response - fitted_values;
                  }
               }
            }
            // Module FIT
            R_vec[tau_index] = Check(residuals_list[tau_index], tau[tau_index]);
            R_lambda_vec[tau_index] = R_vec[tau_index];
            for (k = 0; k < number_penalty_vec[tau_index]; k++)
               R_lambda_vec[tau_index] += lambda * check(sum(jump(_, k) * coefficients_list[tau_index]), tau_penalty);
         }
         // Module CHECK_CONVERGENCE
         R = sum(R_lambda_vec);
         R_lambda = sum(R_lambda_vec);
         if (iter > 10)
            if (std::abs(R_lambda - store_R_lambda) < epsilon_iterations)
               break;
         store_R_lambda = R_lambda;
      }
      for (k = 0; k < number_taus; k++)
      {
         aic_vector[lambda_index] += sample_size * log(R_vec[k]) + 2 * dimension_vec[k];
         bic_vector[lambda_index] += sample_size * log(R_vec[k]) + log(sample_size) * dimension_vec[k];
      }
      ListOf<NumericVector> clone_coefficients_list = clone(coefficients_list);
      ListOf<NumericVector> clone_knots_list = clone(knots_list);
      results[lambda_index] = List::create(_["knots_list"] = clone_knots_list,
                                         _["coefficients_list"] = clone_coefficients_list,
                                         _["dimension_vec"] = dimension_vec);
   }
   results[number_lambdas] = List::create(_["bic_vector"] = bic_vector,
                                          _["aic_vector"] = aic_vector,
                                          _["lambdas"]    = lambdas);
   Rcout << "\n" << "Done!" << std::endl;

   return results;
}
//=========================================================================================
// [[Rcpp::export]]
NumericMatrix bsplines(NumericVector x, NumericVector t,
                       int degree = 1, int derivative = 0)
{
   int x_length = x.size();
   int dimension = t.size() - degree - 1;
   NumericMatrix bm(x_length, dimension);
   for (int j = 0; j < dimension; j++)
      bm(_, j) = bspline(x, t, degree, j, derivative);
   return bm;
}
//=========================================================================================
// [[Rcpp::export]]
NumericVector bspline(NumericVector x, NumericVector t,
                      int degree = 1, int j = 0, int derivative = 0)
{
   int x_length = x.size();
   NumericVector b(x_length);
   int i;
   int k = degree + 1; // k = order = degree + 1
   if (derivative == 0)
   {
      for (i = 0; i < x_length; i++)
      {
         b[i] = 0;
         if ((t[j] <= x[i]) && (x[i] < t[j + k]))
            b[i] = bsp(x[i], t, degree, j);
      }
   }
   else
   {
      for (i = 0; i < x_length; i++)
      {
         b[i] = 0;
         if ((t[j] <= x[i]) && (x[i] < t[j + k]))
            b[i] = dbsp(x[i], t, degree, j, derivative);
      }
   }
   return b;
}
//=========================================================================================
// [[Rcpp::export]]
double bsp(double x, NumericVector t, int degree, int j)
{
   if (degree == 0)
   {
      if ((t[j] <= x) && (x < t[j + 1]))
         return 1;
      else
         return 0;
   }
   else
   {
      double a, b, c, d;
      int k = degree + 1; // k = order = degree + 1
      int jd = j + degree, jpk = j + k, jp1 = j + 1;
      c = t[jd] - t[j];
      if (c > 0)
         a = (x - t[j]) / c;
      else
         a = 0;
      d = t[jpk] - t[jp1];
      if (d > 0)
         b = (t[jpk] - x) / (t[jpk] - t[jp1]);
      else
         b = 0;
      return a * bsp(x, t, degree - 1, j) + b * bsp(x, t, degree - 1, jp1);
   }
}
//=========================================================================================
// [[Rcpp::export]]
double dbsp(double x, NumericVector t, int degree, int j, int derivative)
{
   if (derivative == 0)
      return bsp(x, t, degree, j);
   else
   {
      double a, b, c, d;
      c = t[j + degree] - t[j];
      if (c > 0)
         a = degree / c;
      else
         a = 0;
      int k = degree + 1; // k = order = degree + 1
      int jp1 = j + 1;
      d = t[j + k] - t[jp1];
      if (d > 0)
         b = degree / d;
      else
         b = 0;
      return a * dbsp(x, t, degree - 1, j,     derivative - 1) -
         b * dbsp(x, t, degree - 1, j + 1, derivative - 1);
   }
}
//=========================================================================================
// [[Rcpp::export]]
NumericMatrix jump_bsplines(NumericVector t, int degree = 1)
{
   int k = degree + 1; // k = order = degree + 1
   int dimension = t.size() - k;
   NumericMatrix jump(dimension, dimension - k);
   NumericVector x(dimension - k + 1);
   NumericVector derivative(dimension);
   int j, l;
   for (j = 0; j < x.size(); j++)
      x[j] = 0.5 * (t[j + k - 1] + t[j + k]);
   for (j = 0; j < dimension; j++)
   {
      derivative = bspline(x, t, degree, j, degree);
      for (l = 0; l < (dimension - k); l++)
         jump(j, l) = derivative[l + 1] - derivative[l];
   }
   return jump;
}
//=========================================================================================
// [[Rcpp::export]]
NumericVector knots2t(NumericVector knots, int degree = 1)
{
   double d = mean(diff(knots));
   double min_knots = min(knots);
   double max_knots = max(knots);
   int number_knots = knots.size();
   NumericVector t(2 * degree + number_knots);
   int j;
   for (j = 0; j < degree; j++)
      t[j] = min_knots - d * (degree - j);
   int k = j;
   for (j = 0; j < number_knots; j++)
   {
      t[k] = knots[j];
      k++;
   }
   for (j = 0; j < degree; j++)
   {
      t[k] = max_knots + d * (j + 1);
      k++;
   }
   return t;
}
//=========================================================================================
// [[Rcpp::export]]
NumericVector dim2knots(NumericVector predictor, int dimension, int degree)
{
   int sample_size = predictor.size();
   NumericVector knots(dimension - degree + 1);
   int knot_size = knots.size();
   int j;
   double i;
   for (j = 0; j < (knot_size - 1); j++)
   {
      i = sample_size / (knot_size - 1.0) * j;
      knots[j] = predictor[i];
   }
   knots[knot_size - 1] = max(predictor) + 1e-5;
   return knots;
}
//=========================================================================================
// [[Rcpp::export]]
NumericVector lambdas_all(NumericVector response,
                               int number_lambdas,
                               double lambda_max,
                               double epsilon_lambda)
{
   // compute all lambdas
   // log l_k = log l_1 + (K - k) (log l_K - log l_1)/(K - 1), k = 1, ..., K
   // observe log l_1 = log l_K and log l_K = log l_1
   NumericVector lambdas_all(number_lambdas);
   double lambda_min = epsilon_lambda * lambda_max;
   double ratio_max_min = 1.0 / epsilon_lambda;
   double div, exponent;
   div = number_lambdas - 1;
   for (double lambda_index = 0; lambda_index < number_lambdas; lambda_index++)
   {
      exponent = lambda_index / div;
      lambdas_all[lambda_index] = lambda_min * pow(ratio_max_min, exponent);
   }
   return lambdas_all;
}
//=========================================================================================
// [[Rcpp::export]]
double check(double u, double tau)
{
   if (u < 0)
      return (tau - 1.0) * u;
   else
      return tau * u;
}
//=========================================================================================
// [[Rcpp::export]]
double Check(NumericVector v, double tau)
{
   double risk = 0.0;
   for (int i = 0; i < v.size(); i++)
      risk += check(v[i], tau);
   return risk;
}
//=========================================================================================
// [[Rcpp::export]]
ListOf<NumericVector> NumericVectors(int size)
{
   List lov(size);
   NumericVector v;
   for (int i = 0; i < size; i++)
      lov[i] = v;
   return lov;
}
//=========================================================================================
// [[Rcpp::export]]
ListOf<IntegerVector> IntegerVectors(int size)
{
   List lov(size);
   NumericVector v;
   for (int i = 0; i < size; i++)
      lov[i] = v;
   return lov;
}
//=========================================================================================
// [[Rcpp::export]]
IntegerVector support_of_vector(NumericVector v)
{
   int n = v.size();
   IntegerVector z2n(n);
   for (int i = 0; i < n; i++)
      z2n[i] = i;
   return z2n[abs(v) > 1e-6];
}
//=========================================================================================
// [[Rcpp::export]]
double find_solution(NumericVector v,
                     NumericVector w,
                     double tau,
                     double tau_penalty,
                     double lambda,
                     NumericVector a,
                     NumericVector c,
                     IntegerVector pen_type)
{
   int m = v.size();
   int P = c.size();
   IntegerVector order_c = bubble_order(c);
   a = subsetNumVec(a, order_c);
   pen_type = subsetIntVec(pen_type, order_c);
   c.sort();
   //
   IntegerVector i_star(P);
   for (int p = 0; p < P; p++)
      i_star[p] = find_interval(v, c[p]);
   for (int p = 0; p < P; p++)
   {
      if (!find_slope_v(v, w, tau, tau_penalty, lambda, a, c, i_star[p], pen_type))
         if (find_slope_c(v, w, tau, tau_penalty, lambda, a, c, p, pen_type))
            return c[p];
   }
   int left = 0;
   int right = m - 1;
   int middle;
   while (left < right)
   {
      middle = (right + left) / 2;
      if (!find_slope_v(v, w, tau, tau_penalty, lambda, a, c, middle, pen_type))
         left = middle + 1;
      else
         right = middle;
   }
   return v[left];
}
//=========================================================================================
// [[Rcpp::export]]
bool find_slope_v(NumericVector v,
                  NumericVector w,
                  double tau,
                  double tau_penalty,
                  double lambda,
                  NumericVector a,
                  NumericVector c,
                  int i, IntegerVector pen_type)
{
   int P = c.size();
   double slope;
   double target = v[i];
   if (i == -1)
   {
      target = R_NegInf;
      slope = - tau * sum(w);
   }
   else
      slope = sum(w[seq(0, i)]) - tau * sum(w);
   for (int p = 0; p < P; p++)
   {
      if (pen_type[p] == 1)
         slope += type_forward(target, tau_penalty, lambda, a[p], c[p]);
      else
         slope += type_backward(target, tau_penalty, lambda, a[p], c[p]);
   }
   if (slope >= 0)
      return true;
   else
      return false;
}
//=========================================================================================
// [[Rcpp::export]]
bool find_slope_c(NumericVector v,
                  NumericVector w,
                  double tau,
                  double tau_penalty,
                  double lambda,
                  NumericVector a,
                  NumericVector c,
                  int l, IntegerVector pen_type)
{
   int P = c.size();
   int i_star_l = find_interval(v, c[l]);
   double slope;
   if (i_star_l == -1)
      slope = - tau * sum(w);
   else
      slope = sum(w[seq(0, i_star_l)]) - tau * sum(w);
   for (int p = 0; p < P; p++)
   {
      if (pen_type[p] == 1)
         slope += type_forward(c[l], tau_penalty, lambda, a[p], c[p]);
      else
         slope += type_backward(c[l], tau_penalty, lambda, a[p], c[p]);
   }
   if (slope >= 0)
      return true;
   else
      return false;
}
//=========================================================================================
// [[Rcpp::export]]
double type_forward(double target, double tau_penalty, double lambda, double a, double c)
{
   if (c <= target)
      return lambda * a * tau_penalty;
   else
      return - lambda * a * (1 - tau_penalty);
}
//=========================================================================================
// [[Rcpp::export]]
double type_backward(double target, double tau_penalty, double lambda, double a, double c)
{
   if (c <= target)
      return lambda * a * (1 - tau_penalty);
   else
      return - lambda * a * tau_penalty;
}
//=========================================================================================
// [[Rcpp::export]]
int find_interval(NumericVector y, double z)
{
   int left = 0, right = y.size() - 1;
   while (left <= right)
   {
      int middle = (right + left) / 2;
      // If z greater than y[middle], ignore left half.
      if (y[middle] < z)
         left = middle + 1;
      // If x is smaller, ignore right half.
      else if (y[middle] > z)
         right = middle - 1;
      else
         // z is present at middle
         return middle;
   }
   return (left - 1);
}
//=========================================================================================
// [[Rcpp::export]]
int find_interval_weight(NumericVector w, double tau)
{
   // find k satisfies :
   // w[1:(k-1)] < tau * sum(w) <= w[1:k]
   int left = 0, right = w.size() - 1;
   double tau_w_sum = tau * sum(w);
   // k = 0
   if (tau_w_sum <= w[0])
      return 0;
   // k != 0
   while (left <= right)
   {
      int middle = (right + left) / 2;
      if (sum(w[seq(0, middle)]) < tau_w_sum)
         left = middle + 1;
      else if (sum(w[seq(0, middle)]) > tau_w_sum)
         right = middle - 1;
      else
         return middle;
   }
   return (left);
}
//=========================================================================================
// [[Rcpp::export]]
IntegerVector bubble_order(NumericVector vec)
{
   double tmp = 0;
   NumericVector clone_vec = clone(vec);
   int n = vec.size();
   IntegerVector outvec = seq(1, n);
   int itmp, no_swaps, passes;
   passes = 0;
   while(true)
   {
      no_swaps = 0;
      for (int i = 0; i < n - 1 - passes; ++i)
      {
         if(clone_vec[i] > clone_vec[i+1])
         {
            no_swaps++;
            tmp = clone_vec[i];
            clone_vec[i] = clone_vec[i + 1];
            clone_vec[i + 1] = tmp;
            itmp = outvec[i];
            outvec[i] = outvec[i+1];
            outvec[i+1] = itmp;
         }
      }
      if (no_swaps == 0)
         break;
      passes++;
   }
   return outvec;
}
//=========================================================================================
// [[Rcpp::export]]
NumericVector subsetNumVec(NumericVector x, IntegerVector index)
{
   int n = index.size();
   NumericVector out(n);
   index = index - 1;
   for (int i = 0; i < n; i++)
      out[i] = x[index[i]];
   return out;
}
//=========================================================================================
// [[Rcpp::export]]
IntegerVector subsetIntVec(IntegerVector x, IntegerVector index)
{
   int n = index.size();
   IntegerVector out(n);
   index = index - 1;
   for (int i = 0; i < n; i++)
      out[i] = x[index[i]];
   return out;
}
//=========================================================================================
