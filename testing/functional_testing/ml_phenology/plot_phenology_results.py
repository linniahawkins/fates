#!/usr/bin/env python3
"""
Plot ML Phenology Detection Results

Reads the phenology events and LAI time series from CSV files and generates
visualization plots showing detected Start of Season (SOS) and End of Season (EOS).
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import sys

def plot_phenology_results(output_dir='ml_phenology_output', plot_dir='ml_phenology_plots'):
    """
    Generate plots of LAI time series with detected phenology events.
    
    Parameters:
    -----------
    output_dir : str
        Directory containing the input CSV files
    plot_dir : str
        Directory where plots will be saved
    """
    
    # Create plot directory
    os.makedirs(plot_dir, exist_ok=True)
    
    # Read data
    phenology_file = os.path.join(output_dir, 'phenology_events.csv')
    lai_file = os.path.join(output_dir, 'lai_timeseries.csv')
    
    if not os.path.exists(phenology_file):
        print(f"Error: {phenology_file} not found!")
        sys.exit(1)
    if not os.path.exists(lai_file):
        print(f"Error: {lai_file} not found!")
        sys.exit(1)
    
    phenology = pd.read_csv(phenology_file)
    lai_data = pd.read_csv(lai_file)
    
    print(f"Loaded {len(lai_data)} days of LAI data")
    print(f"Loaded {len(phenology)} years of phenology events")
    
    # Calculate approximate year for each day
    lai_data['year'] = (lai_data['day'] - 1) // 365 + 1
    
    # ==========================================
    # Plot 1: Full time series with all events
    # ==========================================
    fig, ax = plt.subplots(figsize=(16, 6))
    
    ax.plot(lai_data['day'], lai_data['lai'], 'k-', linewidth=0.5, label='LAI')
    
    # Mark SOS and EOS events
    for _, row in phenology.iterrows():
        year = row['year']
        sos_doy = row['sos_doy']
        eos_doy = row['eos_doy']
        
        # Find the day index for this year's events
        year_data = lai_data[lai_data['year'] == year]
        
        if sos_doy > 0:
            sos_days = year_data[year_data['doy'].round(0) == sos_doy]
            if not sos_days.empty:
                day = sos_days.iloc[0]['day']
                lai_val = sos_days.iloc[0]['lai']
                ax.plot(day, lai_val, 'go', markersize=8, 
                       label='SOS' if year == 1 else '', zorder=5)
        
        if eos_doy > 0:
            eos_days = year_data[year_data['doy'].round(0) == eos_doy]
            if not eos_days.empty:
                day = eos_days.iloc[0]['day']
                lai_val = eos_days.iloc[0]['lai']
                ax.plot(day, lai_val, 'ro', markersize=8, 
                       label='EOS' if year == 1 else '', zorder=5)
    
    ax.set_xlabel('Day', fontsize=12)
    ax.set_ylabel('LAI (m²/m²)', fontsize=12)
    ax.set_title('ML-Based Phenology Detection: Full Time Series', fontsize=14, fontweight='bold')
    ax.legend(loc='upper right', fontsize=10)
    ax.grid(True, alpha=0.3)
    
    plot_file = os.path.join(plot_dir, 'phenology_full_timeseries.png')
    plt.tight_layout()
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    print(f"Saved: {plot_file}")
    plt.close()
    
    # ==========================================
    # Plot 2: Individual years (first 4 years)
    # ==========================================
    n_years_to_plot = min(4, len(phenology))
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    axes = axes.flatten()
    
    for i in range(n_years_to_plot):
        ax = axes[i]
        year = i + 1
        year_data = lai_data[lai_data['year'] == year]
        
        if year_data.empty:
            continue
        
        ax.plot(year_data['doy'], year_data['lai'], 'b-', linewidth=1.5, label='LAI')
        
        # Mark SOS
        sos_doy = phenology.loc[phenology['year'] == year, 'sos_doy'].values
        if len(sos_doy) > 0 and sos_doy[0] > 0:
            sos_lai = year_data[year_data['doy'].round(0) == sos_doy[0]]['lai']
            if not sos_lai.empty:
                ax.plot(sos_doy[0], sos_lai.values[0], 'go', markersize=12, 
                       label=f'SOS (DOY {sos_doy[0]})', zorder=5)
                ax.axvline(sos_doy[0], color='g', linestyle='--', alpha=0.5)
        
        # Mark EOS
        eos_doy = phenology.loc[phenology['year'] == year, 'eos_doy'].values
        if len(eos_doy) > 0 and eos_doy[0] > 0:
            eos_lai = year_data[year_data['doy'].round(0) == eos_doy[0]]['lai']
            if not eos_lai.empty:
                ax.plot(eos_doy[0], eos_lai.values[0], 'ro', markersize=12, 
                       label=f'EOS (DOY {eos_doy[0]})', zorder=5)
                ax.axvline(eos_doy[0], color='r', linestyle='--', alpha=0.5)
        
        ax.set_xlabel('Day of Year', fontsize=11)
        ax.set_ylabel('LAI (m²/m²)', fontsize=11)
        ax.set_title(f'Year {year}', fontsize=12, fontweight='bold')
        ax.legend(loc='best', fontsize=9)
        ax.grid(True, alpha=0.3)
        ax.set_xlim(0, 366)
    
    plt.suptitle('ML-Based Phenology Detection: Individual Years', 
                 fontsize=14, fontweight='bold', y=0.995)
    plt.tight_layout()
    plot_file = os.path.join(plot_dir, 'phenology_individual_years.png')
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    print(f"Saved: {plot_file}")
    plt.close()
    
    # ==========================================
    # Plot 3: Phenology statistics
    # ==========================================
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    # Filter out years with no detection (0 values)
    valid_sos = phenology[phenology['sos_doy'] > 0]['sos_doy']
    valid_eos = phenology[phenology['eos_doy'] > 0]['eos_doy']
    
    # SOS distribution
    axes[0].hist(valid_sos, bins=20, color='green', alpha=0.7, edgecolor='black')
    axes[0].axvline(valid_sos.mean(), color='darkgreen', linestyle='--', 
                    linewidth=2, label=f'Mean: {valid_sos.mean():.1f}')
    axes[0].set_xlabel('Day of Year', fontsize=11)
    axes[0].set_ylabel('Frequency', fontsize=11)
    axes[0].set_title('Start of Season (SOS) Distribution', fontsize=12, fontweight='bold')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    
    # EOS distribution
    axes[1].hist(valid_eos, bins=20, color='red', alpha=0.7, edgecolor='black')
    axes[1].axvline(valid_eos.mean(), color='darkred', linestyle='--', 
                    linewidth=2, label=f'Mean: {valid_eos.mean():.1f}')
    axes[1].set_xlabel('Day of Year', fontsize=11)
    axes[1].set_ylabel('Frequency', fontsize=11)
    axes[1].set_title('End of Season (EOS) Distribution', fontsize=12, fontweight='bold')
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plot_file = os.path.join(plot_dir, 'phenology_statistics.png')
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    print(f"Saved: {plot_file}")
    plt.close()
    
    # ==========================================
    # Print summary statistics
    # ==========================================
    print("\n" + "="*50)
    print("PHENOLOGY STATISTICS SUMMARY")
    print("="*50)
    print(f"\nStart of Season (SOS):")
    print(f"  Valid detections: {len(valid_sos)} / {len(phenology)} years")
    if len(valid_sos) > 0:
        print(f"  Mean DOY: {valid_sos.mean():.1f}")
        print(f"  Std Dev: {valid_sos.std():.1f}")
        print(f"  Range: {valid_sos.min():.0f} - {valid_sos.max():.0f}")
    
    print(f"\nEnd of Season (EOS):")
    print(f"  Valid detections: {len(valid_eos)} / {len(phenology)} years")
    if len(valid_eos) > 0:
        print(f"  Mean DOY: {valid_eos.mean():.1f}")
        print(f"  Std Dev: {valid_eos.std():.1f}")
        print(f"  Range: {valid_eos.min():.0f} - {valid_eos.max():.0f}")
    
    # Growing season length
    phenology['gs_length'] = phenology['eos_doy'] - phenology['sos_doy']
    valid_gs = phenology[phenology['gs_length'] > 0]['gs_length']
    if len(valid_gs) > 0:
        print(f"\nGrowing Season Length:")
        print(f"  Mean: {valid_gs.mean():.1f} days")
        print(f"  Std Dev: {valid_gs.std():.1f} days")
        print(f"  Range: {valid_gs.min():.0f} - {valid_gs.max():.0f} days")
    
    print("\n" + "="*50)
    print(f"All plots saved to: {plot_dir}/")
    print("="*50 + "\n")


if __name__ == '__main__':
    # Get directories from command line or use defaults
    output_dir = sys.argv[1] if len(sys.argv) > 1 else 'ml_phenology_output'
    plot_dir = sys.argv[2] if len(sys.argv) > 2 else 'ml_phenology_plots'
    
    print("="*50)
    print("ML Phenology Results Visualization")
    print("="*50)
    print(f"Input directory: {output_dir}")
    print(f"Plot directory: {plot_dir}")
    print("="*50 + "\n")
    
    plot_phenology_results(output_dir, plot_dir)
