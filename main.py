import statistics
from cell import History, M, N, GENE_SIZE, T_MAX


if __name__ == "__main__":
    # Test with smaller values first
    test_mode = False
    if test_mode:
        m, n, t_max = 100, 100, 5
        rates = [0.01]
    else:
        m, n, t_max = M, N, T_MAX
        rates = [0.002, 0.005, 0.01]
    
    filename = "simulation_history"
    for rate in rates:
        print(f"Running simulation with rate: {rate}")
        history = History(m=m, n=n, rate=rate, gene_size=GENE_SIZE, t_max=t_max)
        history.age_multiple_years(t_max)
        # Extract final year data
        final_year_data = history.history[t_max]
        methylation_props = [cell['methylation_proportion'] for cell in final_year_data]
        jsd_scores = [cell['jsd'] for cell in final_year_data]
        
        print(f"Final average methylation proportion for rate {rate}: {statistics.mean(methylation_props):.2%}")
        print(f"Final average JSD for rate {rate}: {statistics.mean(jsd_scores):.4f}")
        history.save_history(filename=filename + f"_{rate:.3f}")