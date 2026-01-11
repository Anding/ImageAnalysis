<#
.SYNOPSIS
    Display histogram visualizations in the terminal using Unicode block characters.

.DESCRIPTION
    This script provides vertical bar chart histogram visualization with sub-character
    resolution for displaying data in the terminal. Each row contains 8 levels of
    resolution using Unicode block characters (▁▂▃▄▅▆▇█).
    
    Can be used in three ways:
    1. Call script directly with parameters: -ScriptValues $data -AutoScale
    2. Dot-source and use Show-Histogram function: . .\Show-Histogram.ps1; Show-Histogram -Values $data
    3. Show examples: .\Show-Histogram.ps1 -Examples

.PARAMETER ScriptValues
    Array of numeric values to display as histogram bars (when calling script directly).

.PARAMETER MaxHeight
    Maximum height in rows for the tallest bar. Default is 20.

.PARAMETER AutoScale
    Automatically scale the data to fit MaxHeight.

.PARAMETER ScaleFactor
    Manual scaling factor if not using AutoScale: BarHeight = Value * ScaleFactor / 8

.PARAMETER Examples
    Display example visualizations.

.EXAMPLE
    # Direct script call with auto-scaling
    .\Show-Histogram.ps1 -ScriptValues @(10, 25, 15, 30, 8) -AutoScale -MaxHeight 10

.EXAMPLE
    # Dot-source and use function
    . .\Show-Histogram.ps1
    Show-Histogram -Values $myHistogram -AutoScale -MaxHeight 15

.EXAMPLE
    # Show built-in examples
    .\Show-Histogram.ps1 -Examples

.EXAMPLE
    # Manual scaling
    .\Show-Histogram.ps1 -ScriptValues @(100, 200, 150) -ScaleFactor 0.5

.NOTES
    Author: GitHub Copilot CLI
    Date: 2026-01-11
    Requires: PowerShell with Unicode support (Windows Terminal recommended)
#>

# Script-level parameters (when called directly)
param(
    [Parameter(Position=0)]
    [double[]]$ScriptValues,
    
    [int]$MaxHeight = 20,
    
    [switch]$AutoScale,
    
    [double]$ScaleFactor = 1.0,
    
    [switch]$Examples
)

function Show-Histogram {
    <#
    .SYNOPSIS
        Displays a vertical bar chart histogram in the terminal using Unicode block characters.
    
    .PARAMETER Values
        Array of numeric values to display as histogram bars.
    
    .PARAMETER MaxHeight
        Maximum height in rows for the tallest bar. Default is 20.
        Each row contains 8 sub-character levels for fine resolution.
    
    .PARAMETER AutoScale
        If specified, automatically scales to MaxHeight based on the maximum value in the data.
        If not specified, uses ScaleFactor parameter.
    
    .PARAMETER ScaleFactor
        Manual scaling factor: BarHeight = Value * ScaleFactor / 8
        Where BarHeight is in rows. Ignored if AutoScale is used.
    
    .EXAMPLE
        Show-Histogram -Values @(10, 25, 15, 30, 8) -AutoScale -MaxHeight 10
    
    .EXAMPLE
        Show-Histogram -Values $histogram -ScaleFactor 0.5 -MaxHeight 20
    #>
    
    [CmdletBinding(DefaultParameterSetName='AutoScale')]
    param(
        [Parameter(Mandatory=$true, Position=0)]
        [double[]]$Values,
        
        [Parameter()]
        [int]$MaxHeight = 20,
        
        [Parameter(ParameterSetName='AutoScale')]
        [switch]$AutoScale,
        
        [Parameter(ParameterSetName='Manual')]
        [double]$ScaleFactor = 1.0
    )
    
    # Block characters for sub-row resolution (0-8 eighths filled)
    $blockChars = " ▁▂▃▄▅▆▇█"
    
    # Calculate heights in eighths (8 levels per row)
    if ($AutoScale) {
        $maxValue = ($Values | Measure-Object -Maximum).Maximum
        if ($maxValue -eq 0) { $maxValue = 1 }  # Avoid division by zero
        $heights = $Values | ForEach-Object { 
            [int](($_ / $maxValue) * $MaxHeight * 8)
        }
    } else {
        $heights = $Values | ForEach-Object { 
            [int]($_ * $ScaleFactor)
        }
    }
    
    $maxHeightInEighths = ($heights | Measure-Object -Maximum).Maximum
    $rowCount = [Math]::Ceiling($maxHeightInEighths / 8)
    
    # Draw from top row to bottom row
    for ($row = $rowCount; $row -gt 0; $row--) {
        $line = ""
        foreach ($height in $heights) {
            $heightAtThisRow = $height - (($row - 1) * 8)
            
            if ($heightAtThisRow -le 0) {
                # Bar doesn't reach this row
                $line += " "
            } elseif ($heightAtThisRow -ge 8) {
                # Full block for this row
                $line += "█"
            } else {
                # Partial block (1-7 eighths)
                $line += $blockChars[$heightAtThisRow]
            }
        }
        Write-Host $line
    }
    
    # Optional: Draw baseline
    Write-Host ("-" * $Values.Count) -ForegroundColor DarkGray
}

# If called directly with parameters, run the function
if ($ScriptValues) {
    if ($AutoScale) {
        Show-Histogram -Values $ScriptValues -AutoScale -MaxHeight $MaxHeight
    } else {
        Show-Histogram -Values $ScriptValues -ScaleFactor $ScaleFactor -MaxHeight $MaxHeight
    }
}
# If called with -Examples or no parameters, show examples
elseif ($Examples -or ($MyInvocation.InvocationName -ne '.' -and -not $ScriptValues)) {
    Write-Host "`nExample 1: Simple test data with AutoScale" -ForegroundColor Cyan
    $testData = @(5, 12, 8, 20, 15, 3, 18, 25, 10, 7)
    Show-Histogram -Values $testData -AutoScale -MaxHeight 10
    
    Write-Host "`nExample 2: Simulated histogram (256 values)" -ForegroundColor Cyan
    # Generate sample histogram data (bell curve-ish)
    $histogram = 0..255 | ForEach-Object {
        [Math]::Exp(-[Math]::Pow(($_ - 128), 2) / 2000) * 1000
    }
    Show-Histogram -Values $histogram -AutoScale -MaxHeight 15
    
    Write-Host "`nExample 3: Manual scaling" -ForegroundColor Cyan
    $moreData = @(100, 250, 180, 320, 150, 90, 280)
    Show-Histogram -Values $moreData -ScaleFactor 0.3 -MaxHeight 20
}
