:: run all estimators, sequentially, for a small number of replicates (and SL library set to SL.glm)
@echo off
setlocal enabledelayedexpansion

:: grab R version for the current machine
For /F "Delims=" %%0 In ('where /r "C:\Program Files\R" Rscript.exe') do set scriptpath="%%~0"
echo Using R executable !scriptpath!

:: pull in command-line argument for scenario
set scenario=%1
set nreps=%2
set estimators=%3
set hal=%4
set link=%5
set ncores=%6
set cgr_type=%7

:: set up arguments to pass
set data_only=0
set use_cached_datasets=0
set S=1
set testing=0
:: send output to directory
set outdir="%cd%\rout"

set found=false
for %%i in (1 3 5 7 9 97 99) do (
  if /i "!scenario!"=="%%i" (
    set found=true
  )
)

IF "!found!"=="true" (
    set "nlist=500 1000 2000 5000"
) ELSE (
    set "nlist=2000 4000 8000 20000"
)

if not exist %outdir% mkdir %outdir%

:: run the simulation
:: loop over n, number of cores
for %%n in (%nlist%) do (
    echo Running scenario %scenario%, estimator %estimators%, n=%%n, cores=%ncores%

  set "this_outfile=%outdir%\output_scenario_%scenario%_link_%link%_nreps%nreps%_n%%n_s%S%_est_%estimators%_cgr_%cgr_type%.out"

  call %scriptpath% 03_main.R ^
      --nreps-total %nreps% ^
      --n %%n ^
      --scenario %scenario% ^
      --estimator %estimators% ^
      --data-only %data_only% ^
      --seed %S% ^
      --use-cached-datasets %use_cached_datasets% ^
      --testing %testing% ^
      --num-cores %ncores% ^
      --eif-reg-type %hal% ^
      --link %link% ^
      --cgr-type %cgr_type% ^
      1>"!this_outfile!" 2>&1
)

