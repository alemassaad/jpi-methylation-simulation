# Statistical Tests Documentation

This document explains the statistical inference tests implemented in the JPI methylation simulation project.

**Implementation**: `phase3/calculate_pvalues.py`

## Overview

The Phase 3 pipeline performs four types of statistical tests to compare methylation metrics across three experimental batches (control, test1, test2):

1. **Student's t-test** - Pairwise comparison assuming equal variance
2. **Welch's t-test** - Pairwise comparison without equal variance assumption
3. **Fisher's ANOVA** - Omnibus test assuming equal variance
4. **Welch's ANOVA** - Omnibus test without equal variance assumption

These tests are applied to 6 comparison metrics:
- Cell JSD Mean
- Cell Methylation Mean
- Gene JSD Mean
- Gene Methylation Mean
- Gene JSD Std Dev
- Gene Methylation Std Dev

---

## Student's t-test

### Explanation
The Student's t-test is a parametric statistical test used to determine whether there is a significant difference between the means of two independent groups. It was developed by William Sealy Gosset under the pseudonym "Student."

### Assumptions
1. **Normality**: Data in each group follows a normal (Gaussian) distribution
2. **Equal Variance (Homoscedasticity)**: The two groups have equal population variances
3. **Independent Samples**: Observations in each group are independent of each other

### Inputs
- **Group 1 data**: Array of numerical values from the first batch (e.g., control)
- **Group 2 data**: Array of numerical values from the second batch (e.g., test1)

In our implementation:
- Each group contains n=3 individuals
- Values represent aggregated metrics (e.g., mean JSD across cells)

### Outputs
- **t-statistic**: Measure of difference between group means relative to variation within groups
- **p-value**: Probability of observing the data (or more extreme) if the null hypothesis (no difference) is true

**Implementation**: `scipy.stats.ttest_ind(data1, data2, equal_var=True)`

---

## Welch's t-test

### Explanation
Welch's t-test is an adaptation of Student's t-test that does not assume equal population variances between the two groups. It is more robust and generally recommended as the default choice for two-sample comparisons.

### Assumptions
1. **Normality**: Data in each group follows a normal distribution
2. **Independent Samples**: Observations in each group are independent of each other
3. **NO equal variance assumption**: Groups can have different variances

### Inputs
- **Group 1 data**: Array of numerical values from the first batch
- **Group 2 data**: Array of numerical values from the second batch

In our implementation:
- Each group contains n=3 individuals
- Values represent aggregated metrics

### Outputs
- **t-statistic**: Modified t-statistic accounting for unequal variances
- **Degrees of freedom**: Adjusted using the Welch-Satterthwaite equation (typically non-integer)
- **p-value**: Probability of observing the data if the null hypothesis is true

**Implementation**: `scipy.stats.ttest_ind(data1, data2, equal_var=False)`

---

## Fisher's ANOVA (One-Way Analysis of Variance)

### Explanation
Fisher's ANOVA, also known as the F-test, is a parametric statistical test used to determine whether there are significant differences among the means of three or more independent groups. It tests the null hypothesis that all group means are equal.

The test compares the variance between groups to the variance within groups. A large F-statistic indicates that between-group differences are larger than expected by chance.

### Assumptions
1. **Normality**: Data in each group follows a normal distribution
2. **Equal Variance (Homoscedasticity)**: All groups have equal population variances
3. **Independent Samples**: Observations within and between groups are independent

### Inputs
- **Group 1 data**: Array of numerical values from control batch (n=3)
- **Group 2 data**: Array of numerical values from test1 batch (n=3)
- **Group 3 data**: Array of numerical values from test2 batch (n=3)

### Outputs
- **F-statistic**: Ratio of between-group variance to within-group variance
- **Degrees of freedom**:
  - **df1** = k - 1 = 2 (number of groups minus one)
  - **df2** = N - k = 6 (total sample size minus number of groups)
- **p-value**: Probability of observing the F-statistic (or larger) if all group means are equal

**Implementation**: `scipy.stats.f_oneway(data_control, data_test1, data_test2)`

**Note**: A significant ANOVA result indicates that at least one group differs from the others, but does not specify which pairs differ. Pairwise t-tests are needed for post-hoc comparisons.

---

## Welch's ANOVA (Heteroscedastic One-Way ANOVA)

### Explanation
Welch's ANOVA is a generalization of Welch's t-test to three or more groups. It tests for differences among group means without assuming equal variances across groups. The test adjusts both the F-statistic calculation and the degrees of freedom to account for heteroscedasticity.

### Assumptions
1. **Normality**: Data in each group follows a normal distribution
2. **Independent Samples**: Observations within and between groups are independent
3. **NO equal variance assumption**: Groups can have different variances

### Inputs
- **Group 1 data**: Array of numerical values from control batch (n=3)
- **Group 2 data**: Array of numerical values from test1 batch (n=3)
- **Group 3 data**: Array of numerical values from test2 batch (n=3)

### Outputs
- **F-statistic**: Modified F-statistic that accounts for unequal variances using weighted group means
- **Degrees of freedom**:
  - **df1** = k - 1 = 2 (number of groups minus one)
  - **df2** = Adjusted using Welch-Satterthwaite approximation (typically non-integer, varies by data)
- **p-value**: Probability of observing the F-statistic (or larger) if all group means are equal

**Implementation**: `scipy.stats.f_oneway(data_control, data_test1, data_test2, equal_var=False)`

**Scipy Version Note**: The `equal_var=False` parameter requires scipy >= 1.15.0. For older versions, the implementation falls back to standard ANOVA.

---

## Pairwise vs Omnibus Tests

### Pairwise Tests (t-tests)
- Compare **two groups at a time**
- Generate 3 comparisons per metric: control vs test1, control vs test2, test1 vs test2
- Total: 3 comparisons × 6 metrics = 18 pairwise tests

### Omnibus Tests (ANOVA)
- Compare **all three groups simultaneously**
- Generate 1 comparison per metric
- Total: 1 comparison × 6 metrics × 2 test types = 12 ANOVA tests

### Multiple Comparisons Correction
With 18 pairwise tests + 12 ANOVA tests = **30 total statistical tests**, the probability of false positives increases. Consider applying Bonferroni correction or other multiple testing corrections when interpreting p-values.

**Bonferroni-corrected significance level**: α_corrected = 0.05 / 30 = 0.00167

---

## Output Files

The statistical tests generate two CSV files in `results/tables/`:

1. **pvalues.csv** (18 rows)
   - Columns: `metric`, `comparison`, `students_t_pvalue`, `welchs_t_pvalue`
   - Contains all pairwise t-test results

2. **anova_results.csv** (12 rows)
   - Columns: `metric`, `test_type`, `f_statistic`, `df1`, `df2`, `p_value`
   - Contains all ANOVA results (Fisher's and Welch's)

---

## Small Sample Size Considerations

Our study design uses **n=3 individuals per batch**, which is a small sample size for parametric tests. Considerations:

1. **Limited power**: Small samples have reduced ability to detect true effects
2. **Normality assumption**: Harder to verify with only 3 observations
3. **Alternative approaches**: Non-parametric tests (Mann-Whitney U, Kruskal-Wallis) may be more appropriate if normality is questionable
4. **Interpretation**: Use p-values cautiously and consider effect sizes alongside significance

Despite these limitations, parametric tests are implemented as they are widely used and provide a foundation for statistical inference in this simulation study.
