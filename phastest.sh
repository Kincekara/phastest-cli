#!/usr/bin/env bash

print_help()
{
    echo "Usage: phastest -i [fasta|genbank|contig] -m [deep|lite] -a [GenBank accession] -s [FASTA sequence file] <option(s)>"
    echo ""
    echo "-i:"
    echo "  Specify input format. Accepted values are 'fasta', 'genbank', or 'contig'."
    echo ""
    echo "-m:"
    echo "  Specify annotation mode. Accepted values are 'deep' or 'lite'."
    echo "  If this flag is not specified, it will default to 'lite'."
    echo "  'deep' uses Prophage Database and PHAST-BSD Bacterial Database."
    echo "  'lite' uses Prophage Database and Swissprot."
    echo "  (Note: 'deep' mode may take significantly longer to complete.)"
    echo ""
    echo "-a:"
    echo "  Specify genbank accession number - parsed only if input-type was set to 'genbank'."
    echo "  e.g.) -a NC_000907.1"
    echo ""
    echo "-s:"
    echo "  Path to the raw FASTA sequence file - parsed only if input-type was set to 'fasta' or 'contig'."
    echo "  Please make sure to provide full FASTA sequence filename that is present within phastest_inputs folder."
    echo "  e.g.) -s {your job}.fna (or .fasta)"
    echo ""
    echo "Options:"
    echo "--yes:"
    echo "  Skip confirmation prompt."
    echo ""
    echo "--silent:"
    echo "  Silence PHASTEST output messages."
    echo ""
    echo "--phage-only:"
    echo "  Only annotate phage region."
}

print_version()
{
    echo "phastest-cli v0.1"
}

## If user sends a termination signal, try to terminate the child job and
## cancel any slurm jobs for the current user, but only if scancel exists.
exit_phastest()
{
    echo "Exiting PHASTEST..."

    # If we started a child job and have its PID, terminate it only.
    if [[ -n "$pid" ]]; then
        if kill -0 "$pid" > /dev/null 2>&1; then
            kill "$pid" 2>/dev/null || true
            sleep 1
            # escalate if still running
            if kill -0 "$pid" > /dev/null 2>&1; then
                kill -9 "$pid" 2>/dev/null || true
            fi
        fi
    fi

    # Only call scancel if it's available on the system.
    if command -v scancel >/dev/null 2>&1; then
        scancel --user="$(whoami)" >/dev/null 2>&1 || true
    fi

    exit 1
}

# Helper function that gets called in for the frontend message.
yes_no()
{
    if [[ $1 == 1 ]]; then
        echo "Yes"
    else
        echo "No"
    fi
}

# Set up the signal trap for common termination signals.
trap exit_phastest SIGINT SIGTERM EXIT

# Set up environments and make sub-programs executable.
# Use PHASTEST_HOME if provided, otherwise default to /phastest-app so
# the installation location can be overridden on different systems.
PHASTEST_HOME=${PHASTEST_HOME:-/phastest-app}

# Only chmod if the directory exists. Avoid assuming absolute paths.
if [[ -d "$PHASTEST_HOME" ]]; then
    chmod -R 755 "$PHASTEST_HOME"
fi

# create JOBS_DIR
cwd=$(pwd)
export JOBS_DIR="$cwd/phastest_out"
if [[ ! -d "$JOBS_DIR" ]]; then
    mkdir "$JOBS_DIR"
    chmod -R 755 "$JOBS_DIR"
fi

input_type=""
anno_mode="lite"
accession=""
sequence=""
filename=""
job_id=""
skip_confirmation=0
silent=0
complete_annotation=1
phage_only=0

# Check for the help flag.
if [[ "$1" == "help" ]]; then
    print_help
    exit 0
fi

# Check for the version flag.
if [[ "$1" == "version" ]]; then
    print_version
    exit 0
fi

# Parse command line arguments.
while getopts ":i:m:a:s:-:" opt; do
    case $opt in
        i)
            input_type=$OPTARG
            ;;
        m)
            anno_mode=$OPTARG
            ;;
        a)
            accession=$OPTARG
            ;;
        s)
            sequence=$OPTARG
            ;;
        -)
            case $OPTARG in
                yes)
                    skip_confirmation=1
                    ;;
                silent)
                    silent=1
                    ;;
                phage-only)
                    complete_annotation=0
                    phage_only=1
                    ;;
                *)
                    echo "Invalid option: --$OPTARG" >&2
                    exit 1
                    ;;
            esac
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            exit 1
            ;;
        :)
            echo "Option -$OPTARG requires an argument." >&2
            exit 1
            ;;
    esac
done

# Check if inputs are all correct.
if [[ $input_type == "genbank" ]] && [[ -z $accession ]]; then
    echo "Genbank accession number is required for 'genbank' input type."
    echo "e.g.) -a NC_000907.1"
    exit 1
fi

if [[ $input_type == "fasta" ]] || [[ $input_type == "contig" ]] && [[ -z $sequence ]]; then
    echo "Sequence file is required for 'fasta' or 'contig' input type."
    echo "e.g.) -s {your job}.fna (or .fasta)"
    exit 1
fi

if [[ $input_type == "" ]]; then
    echo "Input type is required."
    echo "Accepted values are 'fasta', 'genbank', or 'contig'."
    echo "e.g.) -i [fasta|genbank|contig]"
    exit 1
fi

if [[ $anno_mode != "deep" ]] && [[ $anno_mode != "lite" ]]; then
    echo "Invalid annotation mode: $anno_mode"
    echo "Accepted values are 'deep' or 'lite'."
    echo "e.g.) -m [deep|lite]"
    exit 1
fi

# Set up all the parameters for PHASTEST process.
# seq_file_dir="/phastest_inputs"

# if [[ $input_type == "genbank" ]]; then
#     arg1="-g"

#     job_id=$accession
# elif [[ $input_type == "fasta" ]]; then
#     arg1="-s"

#     if [[ ! -f $seq_file_dir/$sequence ]]; then
#         echo "Sequence file not found: $seq_file_dir/$sequence"
#         echo "Please make sure to provide full FASTA sequence filename present within phastest_inputs folder."
#         echo "e.g.) -s seq_test.fna"
#         exit 1
#     fi
# elif [[ $input_type == "contig" ]]; then
#     arg1="-c"

#     if [[ ! -f $seq_file_dir/$sequence ]]; then
#         echo "Sequence file not found: $seq_file_dir/$sequence"
#         echo "Please make sure to provide full FASTA sequence filename present within phastest_inputs folder."
#         echo "e.g.) -s seq_test.fna"
#         exit 1
#     fi
# else 
#     echo "Invalid input type: $input_type"
#     echo "Accepted values are 'fasta', 'genbank', or 'contig'."
#     echo "e.g.) -i [fasta|genbank|contig]"
#     exit 1
# fi

if [[ $input_type == "genbank" ]]; then
    arg1="-g"

    job_id=$accession
elif [[ $input_type == "fasta" ]]; then
    arg1="-s"

    if [[ ! -f "$sequence" ]]; then
        echo "Sequence file not found: $sequence"
        echo "Please make sure to provide full FASTA sequence path."
        #echo "e.g.) -s seq_test.fna"
        exit 1
    fi
elif [[ $input_type == "contig" ]]; then
    arg1="-c"

    if [[ ! -f "$sequence" ]]; then
        echo "Sequence file not found: $sequence"
        echo "Please make sure to provide full FASTA sequence path."
        #echo "e.g.) -s seq_test.fna"
        exit 1
    fi
else 
    echo "Invalid input type: $input_type"
    echo "Accepted values are 'fasta', 'genbank', or 'contig'."
    echo "e.g.) -i [fasta|genbank|contig]"
    exit 1
fi


if [[ $input_type != "genbank" ]]; then
    filename=$(basename "$sequence")
    job_id="${filename%.*}"
fi

# Skip confirmation messages if the "--yes" flag is set.
if [[ $skip_confirmation == 0 ]]; then
    echo "Welcome to PHASTEST!"
    echo ""
    echo "Run PHASTEST with these parameters?"
    echo "Input format: $input_type"
    echo "Job ID: $job_id"
    echo "Annotation Mode: $anno_mode"
    echo ""
    echo "Options:"
    echo "  Skip confirmation prompt (--yes): $(yes_no $skip_confirmation)"
    echo "  Silence PHASTEST output messages (--silent): $(yes_no $silent)"
    echo "  Annotate phage region only (--phage-only): $(yes_no $phage_only)"
    echo ""
    if [[ $sequence != "" ]] && [[ $accession != "" ]]; then
        if [[ $input_type == "genbank" ]]; then
            echo "** Warning: Input type is $input_type; sequence file will be ignored. **"
            echo ""
        else
            echo "** Warning: Input type is $input_type; accession number will be ignored. **"
            echo ""
        fi
    fi
    
    while true; do
        read -p "Continue PHASTEST with these parameters? (Y/N) " yn
        case $yn in
            [Yy]* ) break;;
            [Nn]* ) echo "Exiting PHASTEST..."; exit;;
        esac
    done
fi

if [[ $anno_mode == "deep" ]]; then
    anno_mode="-d"
else
    anno_mode="-l"
fi

# if [[ $input_type == "contig" ]] || [[ $input_type == "fasta" ]]; then
#     mkdir -p "$PHASTEST_HOME/JOBS/$job_id"
#     cp "$sequence" "$PHASTEST_HOME/JOBS/$job_id/$job_id.fna"
# fi

if [[ $input_type == "contig" ]] || [[ $input_type == "fasta" ]]; then
    mkdir -p "$JOBS_DIR/$job_id"
    cp "$sequence" "$JOBS_DIR/$job_id/$job_id.fna"
fi

if [[ $silent == 1 ]]; then
    perl "$PHASTEST_HOME/scripts/phastest.pl" "$arg1" "$job_id" "$anno_mode" "$complete_annotation" > /dev/null 2>&1 &
    pid=$!
else
    perl "$PHASTEST_HOME/scripts/phastest.pl" "$arg1" "$job_id" "$anno_mode" "$complete_annotation" &
    pid=$!
fi

# Wait for the perl script to finish. Receive SIGINT if sent.
wait $pid

exit 0