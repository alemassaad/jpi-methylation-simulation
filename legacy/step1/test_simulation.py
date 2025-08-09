from cell import History

if __name__ == "__main__":
    # Small test with specified parameters
    rate = 0.01   # 1%
    n = 15        # 15 CpG sites per cell
    m = 5         # 5 cells
    t_max = 5     # 5 years
    gene_size = 5 # Default gene size (n=15 is divisible by 5)
    
    print(f"Running small test simulation:")
    print(f"  Rate: {rate:.1%}")
    print(f"  CpG sites per cell (n): {n}")
    print(f"  Number of cells (m): {m}")
    print(f"  Max years (t_max): {t_max}")
    print(f"  Gene size: {gene_size}")
    print()
    
    # Create and run simulation
    history = History(m=m, n=n, rate=rate, gene_size=gene_size, t_max=t_max)
    history.age_multiple_years(t_max)
    
    # Save results
    filename = f"test_small_rate_{rate:.3f}"
    history.save_history(filename, directory="data")
    
    # Display final statistics
    final_year_data = history.history[t_max]
    for i, cell_data in enumerate(final_year_data):
        print(f"Cell {i}: methylation={cell_data['methylation_proportion']:.2%}, JSD={cell_data['jsd']:.4f}")