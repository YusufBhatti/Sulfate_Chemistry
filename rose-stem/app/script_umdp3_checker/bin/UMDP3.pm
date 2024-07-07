# *****************************COPYRIGHT*******************************
# (C) Crown copyright Met Office. All rights reserved. 
# For further details please refer to the file COPYRIGHT.txt 
# which you should have received as part of this distribution. 
# *****************************COPYRIGHT******************************* 

package UMDP3;

# Package to contain subroutines which test for UMDP3 compliance.

# Each subroutine has a standard interface:
# Input:  Array of lines to test
# Output: Scalar value 0=pass, >0 = fail

# Subroutines which don't obey this interface:
#  get_include_number - returns the value of a variable scoped to this file
#                   (number of files using includes for variable declarations)
#  remove_quoted      - returns the input string having removed any quoted
#                        substrings (single or double).

# Standard modules
use strict;
use warnings;
use 5.010;

use Text::Balanced qw(extract_quotelike extract_multiple);

# Global variables 

my $number_of_files_with_variable_declarations_in_includes = 0;

sub get_include_number {
  return $number_of_files_with_variable_declarations_in_includes;
}

sub remove_quoted {
  my $line = shift;

  # Replace quoted strings with a blessed reference:
  my @strings = extract_multiple($line, [
                                 {Quoted =>  sub { extract_quotelike($_[0]) } },
                                 ] );
  # Stitch the non-quoted fields back together into a single string:
  my $remainder = "";
  foreach my $string (@strings) {
    $remainder .= $string if not ($string =~ /^Quoted=SCALAR/);
  }
  return $remainder;
}


my @fortran_keywords = (
    '\.FALSE\.',
    '\.TRUE\.',
    'ABORT',
    'ABS',
    'ACCESS',
    'ACHAR',
    'ACOS',
    'ACOSH',
    'ADJUSTL',
    'ADJUSTR',
    'AIMAG',
    'AINT',
    'ALARM',
    'ALGAMA',
    'ALL',
    'ALLOCATABLE',
    'ALLOCATED',
    'ALOG',
    'ALOG10',
    'AMAX0',
    'AMAX1',
    'AMIN0',
    'AMIN1',
    'AMOD',
    '\.AND\.',
    'ANINT',
    'ANY',
    'ASIN',
    'ASINH',
    'ASSOCIATED',
    'ATAN',
    'ATAN2',
    'ATANH',
    'ATOMIC_DEFINE',
    'ATOMIC_REF',
    'BESJ0',
    'BESJ1',
    'BESJN',
    'BESSEL_J0',
    'BESSEL_J1',
    'BESSEL_JN',
    'BESSEL_Y0',
    'BESSEL_Y1',
    'BESSEL_YN',
    'BESY0',
    'BESY1',
    'BESYN',
    'BGE',
    'BGT',
    'BIT_SIZE',
    'BLE',
    'BLT',
    'BTEST',
    'CABS',
    'CALL',
    'CCOS',
    'CDABS',
    'CDCOS',
    'CDEXP',
    'CDLOG',
    'CDSIN',
    'CDSQRT',
    'CEILING',
    'CEXP',
    'CHAR',
    'CHARACTER',
    'CHDIR',
    'CHMOD',
    'CLOG',
    'CMPLX',
    'COMMAND_ARGUMENT_COUNT',
    'COMPILER_OPTIONS',
    'COMPILER_VERSION',
    'COMPLEX',
    'CONJG',
    'CONTAINS',
    'CONVERT',
    'COS',
    'COSH',
    'COUNT',
    'CPP',
    'CPU_TIME',
    'CSHIFT',
    'CSIN',
    'CSQRT',
    'CTIME',
    'CYCLE',
    'C_ASSOCIATED',
    'C_FUNLOC',
    'C_F_POINTER',
    'C_F_PROCPOINTER',
    'C_LOC',
    'C_SIZEOF',
    'DABS',
    'DACOS',
    'DACOSH',
    'DASIN',
    'DASINH',
    'DATAN',
    'DATAN2',
    'DATANH',
    'DATE_AND_TIME',
    'DBESJ0',
    'DBESJ1',
    'DBESJN',
    'DBESY0',
    'DBESY1',
    'DBESYN',
    'DBLE',
    'DCMPLX',
    'DCONJG',
    'DCOS',
    'DCOSH',
    'DDIM',
    'DECODE',
    'DEXP',
    'DFLOAT',
    'DGAMMA',
    'DIGITS',
    'DIM',
    'DIMAG',
    'DINT',
    'DLGAMA',
    'DLOG',
    'DLOG10',
    'DMAX1',
    'DMIN1',
    'DMOD',
    'DNINT',
    'DO',
    'DOT_PRODUCT',
    'DPROD',
    'DREAL',
    'DSHIFTL',
    'DSHIFTR',
    'DSIGN',
    'DSIN',
    'DSINH',
    'DSQRT',
    'DTAN',
    'DTANH',
    'DTIME',
    'ENCODE',
    'END DO',
    'END IF',
    'ENUM',
    'ENUMERATOR',
    'EOSHIFT',
    'ELSE',
    'EPSILON',
    'ERF',
    'ERFC',
    'ERFC_SCALED',
    'ETIME',
    'EXECUTE_COMMAND_LINE',
    'EXIT',
    'EXP',
    'EXPONENT',
    'EXTENDS_TYPE_OF',
    'FDATE',
    'FGET',
    'FGETC',
    'FLOAT',
    'FLOOR',
    'FLUSH',
    'FNUM',
    'FORALL',
    'FORMAT',
    'FPP',
    'FPUT',
    'FPUTC',
    'FRACTION',
    'FREE',
    'FSEEK',
    'FSTAT',
    'FTELL',
    'FUNCTION',
    'GAMMA',
    'GERROR',
    'GETARG',
    'GETCWD',
    'GETENV',
    'GETGID',
    'GETLOG',
    'GETPID',
    'GETUID',
    'GET_COMMAND',
    'GET_COMMAND_ARGUMENT',
    'GET_ENVIRONMENT_VARIABLE',
    'GMTIME',
    'HOSTNM',
    'HUGE',
    'HYPOT',
    'IABS',
    'IACHAR',
    'IALL',
    'IAND',
    'IANY',
    'IARGC',
    'IBCLR',
    'IBITS',
    'IBSET',
    'ICHAR',
    'IDATE',
    'IDIM',
    'IDINT',
    'IDNINT',
    'IEOR',
    'IERRNO',
    'IF',
    'IFIX',
    'IMAG',
    'IMAGE_INDEX',
    'IMAGPART',
    'IMPLICIT',
    'IMPORT',
    'IN',
    'INCLUDE',
    'INDEX',
    'INT',
    'INT2',
    'INT8',
    'INTEGER',
    'INTENT',
    'INTERFACE',
    'INTRINSIC',
    'IOR',
    'IPARITY',
    'IRAND',
    'ISATTY',
    'ISHFT',
    'ISHFTC',
    'ISIGN',
    'ISNAN',
    'ISO_FORTRAN_ENV',
    'IS_IOSTAT_END',
    'IS_IOSTAT_EOR',
    'KILL',
    'KIND',
    'LBOUND',
    'LCOBOUND',
    'LEADZ',
    'LEN',
    'LEN_TRIM',
    'LGAMMA',
    'LGE',
    'LGT',
    'LINK',
    'LLE',
    'LLT',
    'LNBLNK',
    'LOC',
    'LOG',
    'LOG10',
    'LOGICAL',
    'LOG_GAMMA',
    'LONG',
    'LSHIFT',
    'LSTAT',
    'LTIME',
    'MALLOC',
    'MASKL',
    'MASKR',
    'MATMUL',
    'MAX',
    'MAX0',
    'MAX1',
    'MAXEXPONENT',
    'MAXLOC',
    'MAXVAL',
    'MCLOCK',
    'MCLOCK8',
    'MERGE',
    'MERGE_BITS',
    'MIN',
    'MIN0',
    'MIN1',
    'MINEXPONENT',
    'MINLOC',
    'MINVAL',
    'MOD',
    'MODULE',
    'MODULO',
    'MOVE_ALLOC',
    'MVBITS',
    'NAMELIST',
    'NEAREST',
    'NEW_LINE',
    'NINT',
    'NONE',
    'NORM2',
    '\.NOT\.',
    'NULL',
    'NULLIFY',
    'NUM_IMAGES',
    'OPERATOR',
    'OPTIONAL',
    '\.OR\.',
    'OUT',
    'PACK',
    'PARITY',
    'PERROR',
    'POINTER',
    'POPCNT',
    'POPPAR',
    'PRECISION',
    'PRESENT',
    'PRIVATE',
    'PRODUCT',
    'PROTECTED',
    'PUBLIC',
    'RADIX',
    'RAN',
    'RAND',
    'RANDOM_NUMBER',
    'RANDOM_SEED',
    'RANGE',
    'RANK',
    'READ',
    'REAL',
    'REALPART',
    'RECORD',
    'RENAME',
    'REPEAT',
    'RESHAPE',
    'RRSPACING',
    'RSHIFT',
    'SAME_TYPE_AS',
    'SAVE',
    'SCALE',
    'SCAN',
    'SECNDS',
    'SECOND',
    'SELECTED_CHAR_KIND',
    'SELECTED_INT_KIND',
    'SELECTED_REAL_KIND',
    'SET_EXPONENT',
    'SHAPE',
    'SHIFTA',
    'SHIFTL',
    'SHIFTR',
    'SHORT',
    'SIGN',
    'SIGNAL',
    'SIN',
    'SINH',
    'SIZE',
    'SIZEOF',
    'SLEEP',
    'SNGL',
    'SPACING',
    'SPREAD',
    'SQRT',
    'SRAND',
    'STAT',
    'STORAGE_SIZE',
    'STRUCTURE',
    'SUBROUTINE',
    'SUM',
    'SYMLNK',
    'SYSTEM',
    'SYSTEM_CLOCK',
    'TAN',
    'TANH',
    'THIS_IMAGE',
    'THEN',
    'TIME',
    'TIME8',
    'TINY',
    'TRAILZ',
    'TRANSFER',
    'TRANSPOSE',
    'TRIM',
    'TTYNAM',
    'TYPE',
    'UBOUND',
    'UCOBOUND',
    'UMASK',
    'UNLINK',
    'UNPACK',
    'USE',
    'VALUE',
    'VERIFY',
    'VOLATILE',
    'WHERE',
    'WRITE',
    '\.XOR\.',
    'ZABS',
    'ZCOS',
    'ZEXP',
    'ZLOG',
    'ZSIN',
    'ZSQRT',
);


my @archaic_fortran_keywords = (
    'ALOG',
    'ALOG10',
    'AMAX0',
    'AMAX1',
    'AMIN0',
    'AMIN1',
    'AMOD',
    'CABS',
    'CCOS',
    'CEXP',
    'CLOG',
    'CSIN',
    'CSQRT',
    'DABS',
    'DACOS',
    'DASIN',
    'DATAN',
    'DATAN2',
    'DBESJ0',
    'DBESJ1',
    'DBESJN',
    'DBESY0',
    'DBESY1',
    'DBESYN',
    'DCOS',
    'DCOSH',
    'DDIM',
    'DERF',
    'DERFC',
    'DEXP',
    'DINT',
    'DLOG',
    'DLOG10',
    'DMAX1',
    'DMIN1',
    'DMOD',
    'DNINT',
    'DSIGN',
    'DSIN',
    'DSINH',
    'DSQRT',
    'DTAN',
    'DTANH',
    'FLOAT',
    'IABS',
    'IDIM',
    'IDINT',
    'IDNINT',
    'IFIX',
    'ISIGN',
    'LONG',
    'MAX0',
    'MAX1',
    'MIN0',
    'MIN1',
    'SNGL',
    'ZABS',
    'ZCOS',
    'ZEXP',
    'ZLOG',
    'ZSIN',
    'ZSQRT',
    );
                       
my @openmp_keywords = (
    'PARALLEL', 'MASTER','CRITICAL', 'ATOMIC', 'SECTIONS', 'WORKSHARE', 'TASK',
    'BARRIER', 'TASKWAIT', 'FLUSH', 'ORDERED', 'THREADPRIVATE', 'SHARED', 
    'DEFAULT', 'FIRSTPRIVATE', 'LASTPRIVATE', 'COPYIN', 'COPYPRIVATE', 
    'REDUCTION',
                       );

my @fortran_types = (
    'TYPE',
    'CLASS',
    'INTEGER',
    'REAL',
    'DOUBLE PRECISION',
    'CHARACTER',
    'LOGICAL',
    'COMPLEX',
    'ENUMERATOR',
                       );

sub get_fortran_keywords {
  return @fortran_keywords;
}

sub get_openmp_keywords {
  return @openmp_keywords;
}

sub get_archaic_fortran_keywords {
  return @archaic_fortran_keywords;
}

# List of uncapitalised keywords present in most recently tested file
my %extra_error_information = ();
sub get_extra_error_information {
  return %extra_error_information;
}

sub reset_extra_error_information {
  %extra_error_information = ();
}
################################# UMDP3 tests #################################


# Check for uncapitalised keywords
sub capitalised_keywords {
  my @lines = @_;
  
 
  my $failed = 0;

# Iterate over lines and keywords                 
  foreach my $line (@lines) {
    my @keywords_to_check = get_fortran_keywords();

    $line = remove_quoted($line);
    next unless $line;
    next unless $line =~/\S/;   # If line empty, try the next

# Remove comments unless they're OpenMP commands or unless they're in quotes
    if ($line =~/![^\$]/) {
      unless ($line =~/".*!.*"/ or $line =~/'.*!.*'/) {
        $line =~s/![^\$].*//g;
      }
    }
    
    if ($line =~/^!\$/) {
      push @keywords_to_check, get_openmp_keywords();
    }
    
    foreach my $keyword (@keywords_to_check) {


# If the keyword is present on the line
      if ($line =~/\W$keyword\W/i or $line =~/^$keyword\W/i or 
          $line =~/\%$keyword\W /i or $line =~/\($keyword\)/i) { 
# If the line contains WRITE let it pass as it's probably text
        if ($line =~/WRITE.*'/ or $line =~/WRITE.*"/ ) {
          next;
        }

        if ($line =~/\(kind\s*=.*::/) {
          $extra_error_information{'KIND'}++;
          $failed++;
        }
# If the keyword is followed by an underscore it's a variable, not a keyword
        if ($line =~/$keyword\_/i) {
          next;
        }

        # Ignore keyword between quotes
        next if ($line =~/'.*$keyword.*'/i);
        next if ($line =~/".*$keyword.*"/i);

        # Ignore cases such as RESHAPE(len=something) where 'len' would
        # otherwise be triggered
        next if ($line =~/,$keyword=/i);
        next if ($line =~/\($keyword=/i);

# If the keyword is not uppercase  
        unless ($line =~ /$keyword\W/) {
# Ignore CPP
          unless ($line =~/^\s*#/) {
            $extra_error_information{$keyword}++;
            $failed++;
          }
        }
      }
    }
  }

  return $failed;
}


# OpenMP sentinels must be in column one
sub openmp_sentinels_in_column_one {
  my @lines = @_;
  
  my $failed = 0;
  foreach my $line (@lines) {
# Check for one or more spaces before !$
    $failed++ if ($line =~/\s+!\$/);
  }

  return $failed;
}

# ENDIF, etc should be END IF
sub unseparated_keywords {
  my @lines = @_;


  my @unseparated_keywords = (
    'ENDIF', 'ENDDO', 'ENDWHERE', 'ELSEIF', 'ENDINTERFACE',
    'ENDSUBROUTINE', 'ENDMODULE', 'ENDFUNCTION', 'ENDSELECT', 'ENDPARALLEL',
    'ENDPARALLELDO', 'ENDCRITICAL', 'PARALLELDO',
                             );
    
  
  my $failed = 0;
  foreach my $line (@lines) {

# Remove comments unless they're OpenMP commands
    if ($line =~/![^\$]/) {
        $line =~s/![^\$].*//g;
    }

# Check for frequent ones - should rewrite as a loop
    unless ($line =~/^\s*#/) {     # Ignore CPP
      foreach my $keyword (@unseparated_keywords) {
        if ($line =~/\W$keyword\W/ or $line=~/^$keyword/) {
          $failed++;
          $extra_error_information{$keyword}++;
        }
      }
    }
  }

  return $failed;
}

# PAUSE and EQUIVALENCE are forbidden
sub forbidden_keywords {
  my @lines = @_;
  
  my $failed = 0;
  foreach my $line (@lines) {

# Remove comments unless they're OpenMP commands
    if ($line =~/![^\$]/) {
        $line =~s/![^\$].*//g;
    }

    $failed++ if ($line =~/\WEQUIVALENCE\W/i);
    $failed++ if ($line =~/\WPAUSE\W/i);
    $failed++ if ($line =~/^EQUIVALENCE\W/i);
    $failed++ if ($line =~/^PAUSE\W/i);
  }

  return $failed;
}

# Older forms of relational operators are forbidden
sub forbidden_operators {
  my @lines = @_;
  
  my $failed = 0;
  foreach my $line (@lines) {

# Remove comments unless they're OpenMP commands
    if ($line =~/![^\$]/) {
        $line =~s/![^\$].*//g;
    }

    $failed++ if ($line =~/\.GT\./i);
    $failed++ if ($line =~/\.GE\./i);
    $failed++ if ($line =~/\.LT\./i);
    $failed++ if ($line =~/\.LE\./i);
    $failed++ if ($line =~/\.EQ\./i);
    $failed++ if ($line =~/\.NE\./i);
  }

  return $failed;
}


# Any GO TO must go to 9999
sub go_to_other_than_9999 {
  my @lines = @_;
  
  my $failed = 0;
  foreach my $line (@lines) {

# Remove comments unless they're OpenMP commands
    if ($line =~/![^\$]/) {
        $line =~s/![^\$].*//g;
    }
# Delete quoted strings
    $line =~s/".*"//;
    $line =~s/'.*'//;
    
# Find lines matching GO TO
    if ($line =~/GO\s*TO/i) {
# If the line number isn't 9999
      unless ($line =~/GO\s*TO\s*9999/i) {
        $failed++;
      }
    } 

  }

  return $failed;
}

# WRITE must specify a proper format
sub write_using_default_format {
  my @lines = @_;
  
  my $failed = 0;
  foreach my $line (@lines) {

# Remove comments unless they're OpenMP commands
    if ($line =~/![^\$]/) {
        $line =~s/![^\$].*//g;
    }

# Check for WRITE(...*)
    if ($line =~/WRITE\s*\(.*\*\)/i) {
      $failed++;
    } 

  }

  return $failed;
}

sub lowercase_variable_names {
  my @lines = @_;
  
  my $failed = 0;
  my @variables;


# Make a list of variables
  foreach my $line (@lines) {

# Remove comments unless they're OpenMP commands
    if ($line =~/![^\$]/) {
        $line =~s/![^\$].*//g;
    }

    if ($line =~/^\s*REAL/i or $line =~/^\s*INTEGER/i or $line =~/^\s*LOGICAL/ or $line=~/^\s*CHARACTER/) {
      if ($line =~/::/) {
        $line =~/::\s*(\w+)/;
        my $variable = $1;
        next unless ($variable);
        push @variables, $variable;
      }
    }
  }
  
# Search the code for these variables
  foreach my $line (@lines) {
    # Ignore CPP defs:
    next if ($line =~/^\s*#/);
    foreach my $variable (@variables) {
      if ($line =~/\b($variable)\b/i) {
        my $instance_of_variable = $1;
        # If the variable is 4 or more characters and is uppercase in the declaration fail the test 
        # The length test is because some short scientific quantities could legitimately be uppercase.
        next if (length $variable < 4 );
        
        # Ignore variables in quotes
        next if ($line =~/'.*$variable.*'/i);
        next if ($line =~/".*$variable.*"/i);
        
        if ($instance_of_variable eq "\U$instance_of_variable") {
          $failed++;
          $extra_error_information{$instance_of_variable}++;
        }
      }
    }
  }
  
  return $failed;  
}

sub include_files_for_variable_declarations {
  my @lines = @_;
  
  my $failed = 0;

  my $found_dr_hook = 0;
  foreach my $line (@lines) {
    $found_dr_hook++ if ($line =~/CALL\s+dr_hook/i);
  }

  # File which don't have directly executable code automatically pass this
  return 0 unless $found_dr_hook;

  foreach my $line (@lines) {
    $failed++ if ($line =~/^\s*#include/);
    last if ($line=~/CALL\s+dr_hook/i);
  }

  $number_of_files_with_variable_declarations_in_includes++ if $failed;  
  return $failed;
}

sub dimension_forbidden {
  my @lines = @_;
  
  my $failed = 0;
  foreach my $line (@lines) {

# Remove comments unless they're OpenMP commands
    if ($line =~/![^\$]/) {
        $line =~s/![^\$].*//g;
    }
    $line = remove_quoted($line);
    next unless $line;
    next if ($line =~/WRITE.*DIMENSION/i);
    next if ($line =~/_DIMENSION/i);
    next if ($line =~/DIMENSION_/i);
    $failed++ if ($line =~/\WDIMENSION\W/i);
    $failed++ if ($line =~/^DIMENSION\W/i);
  }
  
  return $failed;
}

sub forbidden_stop {
  my @lines = @_;
  
  my $failed = 0;
  foreach my $line (@lines) {

# Remove comments unless they're OpenMP commands
    if ($line =~/![^\$]/) {
        $line =~s/![^\$].*//g;
    }

    $failed++ if ($line =~/^\s*STOP\s/i);
    $failed++ if ($line =~/^\s*CALL\s*abort\W/i);
  }
  
  return $failed;
}

sub ampersand_continuation {
  my @lines = @_;
  
  my $failed = 0;
  foreach my $line (@lines) {

# Remove comments unless they're OpenMP commands
    if ($line =~/![^\$]/) {
        $line =~s/![^\$].*//g;
    }

    $failed++ if ($line =~/^\s*&/i);
    $failed++ if ($line =~/^\s*!\$\s*&/i);
  }
  
  return $failed;
}

sub implicit_none {
  my @lines = @_;
  
  my $failed = 0;
  my $foundit = 0;
  my $modules = 0;
  my @lines_to_test;
  
  my $in_interface = 0;
  foreach my $input_line (@lines) {
# Remove comments unless they're OpenMP commands
    if ($input_line =~/![^\$]/) {
        $input_line =~s/![^\$].*//g;
    }
# MODULEs etc in INTERFACEs don't have implicit none, so ignore these
    if ($input_line =~/^\s*INTERFACE\s/i) {
      $in_interface = 1;
    }
    push @lines_to_test, $input_line unless $in_interface;
    if ($input_line =~/^\s*END\s*INTERFACE/i) {
      $in_interface = 0;
    }
  }
  
  foreach my $line (@lines_to_test) {

    $foundit++ if ($line =~/^\s*IMPLICIT\s+NONE/i);
    $modules++ if ($line =~/^\s*SUBROUTINE\W/i or $line =~/^\s*MODULE\W/i or 
                   $line =~/^\s*FUNCTION\W/i or 
                   $line =~/^\s*REAL\s*FUNCTION\W/i or
                   $line =~/^\s*LOGICAL\s*FUNCTION\W/i or 
                   $line =~/^\s*INTEGER\s*FUNCTION\W/i or 
                   $line =~/^\s*PROGRAM\W/i);
  }
  
  $failed = 1 unless ($foundit >= $modules);
  
  return $failed;
}

sub intrinsic_as_variable {
  my @lines = @_;
  my $failed = 0;
  my @keywords = get_fortran_keywords();

  my @fixed_lines=();

  # Steps:
  #  i)   sanitise lines
  #  ii)  look for match
  #  iii) check if match is a declaration (which must start with a type)
  #  iv)  exclude any false positives from initialisation.

  # i) sanitise lines

  foreach my $line (@lines) {
    # Remove comments unless they're OpenMP commands or unless they're in quotes
    if ($line =~/![^\$]/) {
      unless ($line =~/"[^"]*?![^"]*?"/ or $line =~/'[^']*?![^']*?'/) {
        $line =~s/![^\$].*//g;
      }
    }
    $line = remove_quoted($line);

    # Remove pre-processing directives
    if ($line =~/\s*#/) {
      $line = "";
    }

    push @fixed_lines, $line;
  }

  my $entire = join("",@fixed_lines);

  # Sort out continuation lines
  $entire =~s/&\s*\n//g;

  @fixed_lines = split/\n/,$entire;

  foreach my $line (@fixed_lines) {

    next unless $line;
    next unless $line =~/\S/;

    my $oline = $line;

    foreach my $keyword (@keywords) {
      my $decl_match = 0;
      $line = $oline;

      #  ii)  look for match
      if ($line =~/(^|\W)$keyword($|\W)/i) {
        foreach my $type (@fortran_types) {

          #  iii) check if match is a variable declaration (which always starts with a type):
          if ($line =~/^\s*$type(\W.*\W|\W)$keyword/i) {
              $decl_match++;
              last;
          }
        }
      }

      #  iv)  exclude any false positives from initialisation.
      if ($decl_match > 0) {
        # This is a variable declaration with a matching keyword
        # make sure this is not because of initialising to the result of a function
        # (i.e. the keyword is the RHS of the = in this variable initialisation).

        # remove any type attributes which may match the keyword
        $line =~s/^.*:://g;

        # at this point, things in brackets aren't relevant
        while ($line =~ /\(.*\)/) {
          $line =~s/\([^()]*?\)//g;
        }

        # split on commas, in case there are multiple variable declarations
        my @decls = split/,/,$line;

        foreach my $decl (@decls) {

          # RHS of = signs are not variable definition (are instead initialiser etc.)
          $decl =~s/=.*$//g;

          # Remove function declarations
          $decl =~s/^.*?\sFUNCTION\s//ig;

          # If we get this far any matches are fails
          if ($decl =~/(^|\W)$keyword(\W|$)/i) {
            $line = "\n    $keyword";
            $failed++;
            $extra_error_information{$line}++;
          }
        }
      }
    }
  }

  return $failed;
}

sub line_over_80chars {
  my @lines = @_;
  my $failed = 0;

  foreach my $line (@lines) {
    # This needs to be 81, as Perl counts the newline as having length 1
    if (length $line > 81) {
      $failed++;
      # Reformat line so it prints the offending line neatly
      chomp($line);
      $line = "'$line'";
      $extra_error_information{$line}++;
    }
  }
  return $failed;
}

sub tab_detection {
  my @lines = @_;
  my $failed = 0;
  foreach my $line (@lines) {
    # If any line contains a tab character
    $failed++ if ($line=~/\t/);
  }
  return $failed;
}

sub check_crown_copyright {
  my @lines = @_;
  my $failed = 1;
  my @valid_agreements = ( 
        'L0195', 'NERC',      'SC0138', 'UKCA',
        'SC0171', 'ACCESS',   'SC0237', 'JULES',
        'IBM',
                        );
  
  foreach my $line (@lines) {
    $failed = 0 if ($line =~/^\s*(!|\/\*).*Crown\s*copyright/i);
    foreach my $agreement (@valid_agreements) {
      $failed = 0 if ($line =~/^\s*(!|\/\*).*$agreement/i);
    }
  }
  
  return $failed;
}


sub check_code_owner {
  my @lines = @_;
  my $failed = 1;
  my $failed_co = 0;
  my $failed_bi = 0;

  foreach my $line (@lines) {
      $failed_co++ if ($line=~ /^\s*(!|\/\*)\s*Code\s*Owner:\s*Please\s*refer\s*to\s*the\s*UM\s*file\s*CodeOwners\.txt/i);
      $failed_bi++ if ($line=~ /^\s*(!|\/\*)\s*This\s*file\s*belongs\s*in\s*section:/i);
    }

    if ($failed_co > 1 or $failed_bi > 1) {
      $extra_error_information{"(multiple statements found)"}++;
    }


    if ($failed_co == 1 and $failed_bi == 1) {
      $failed=0;
    }

  return $failed;
}


sub retire_if_def {
  my @lines = @_;
  my @ifdefs = (
    'VATPOLES', 'A12_4A', 'A12_3A', 'UM_JULES', 'A12_2A',
               );

  # Sort out C continuation lines
  my $entire = join("",@lines);
  $entire =~s/\\\s*\n//g;
  @lines = split/\n/,$entire;

  my $failed = 0;
  for (my $i=0; $i < scalar @lines; $i++) {
    my $line = $lines[$i];
    foreach my $ifdef (@ifdefs) {
      # matches #if defined(<def>), #elif defined(<def>), #ifdef <def>, and #ifndef <def>
      if ($line =~/^\s*#(el)?if.*\W$ifdef/) {
        $failed++;
        $extra_error_information{$ifdef}++;
      }
    }
  }

  return $failed;
}

sub c_deprecated {
  my @lines = @_;
  my %deprecateds = (
    'strcpy' => "(): please use strncpy() instead",
    'sprintf' => "(): please use snprintf() instead",
    );

  my $entire = join("",@lines);

  #remove commented sections
  $entire =~s/\/\*(.|\n)+?(\*\/)//g;

  # Sort out continuation lines
  $entire =~s/\\\s*\n//g;

  #remove #pragmas
  $entire =~s/(^|\n)\s*#pragma.+?\n/\n/g;

  @lines = split/\n/,$entire;

  my $failed = 0;
  for (my $i=0; $i < scalar @lines; $i++) {
    my $line = $lines[$i];

    foreach my $dep (keys %deprecateds) {
      if ($line =~/$dep/) {
        my $extra_msg = "$dep$deprecateds{$dep}" ;
        $failed++;
        $extra_error_information{$extra_msg}++;
      }
    }
  }
    
  return $failed;
}

sub printstatus_mod {
  my @lines = @_;
  
  my $failed = 0;
  foreach my $line (@lines) {
    $failed++ if ($line =~/^\s*USE\s*printstatus_mod/i);
  }
  return $failed;
}

sub write6 {
  my @lines = @_;
  my $failed = 0;
  
  for (my $i = 0; $i < scalar @lines; $i++) {
    if ($lines[$i] =~ /^\s*WRITE/i) {
      if ($lines[$i] =~/^\s*WRITE\s*\(\s*6/) {
        $failed++;
      }
    }
    
  }

  return $failed;
}

sub printstar {
  my @lines = @_;
  
  my $failed = 0;
  foreach my $line (@lines) {
    $failed++ if ($line =~/^\s*PRINT\s*\*/i);
  }
  return $failed;
}
  
sub um_fort_flush {
  my @lines = @_;
  
  my $failed = 0;
  foreach my $line (@lines) {
    $failed++ if ($line =~/^\s*CALL\s*UM_FORT_FLUSH/i);
  }
  return $failed;
}

sub svn_keyword_subst {
  my @lines = @_;
 
  my $failed = 0;
  foreach my $line (@lines) {
    $failed++ if ($line =~/\$Date\$/);
    $failed++ if ($line =~/\$LastChangedDate\$/);
    $failed++ if ($line =~/\$Revision\$/);
    $failed++ if ($line =~/\$Rev\$/);
    $failed++ if ($line =~/\$LastChangedRevision\$/);
    $failed++ if ($line =~/\$Author\$/);
    $failed++ if ($line =~/\$LastChangedBy\$/);
    $failed++ if ($line =~/\$HeadURL\$/);
    $failed++ if ($line =~/\$URL\$/);
    $failed++ if ($line =~/\$Id\$/);
    $failed++ if ($line =~/\$Header\$/);
  }
  
  return $failed;

}

sub omp_missing_dollar {
  my @lines = @_;
 
  my $failed = 0;
  foreach my $line (@lines) {
    if ($line =~/^\s*!OMP/) {
      $failed = 1;
    }
  }
  
  return $failed;

}

sub cpp_ifdef {
  # ifdefs should be of the form "#if defined(MY_IFDEF)"
  # rather than "#ifdef(MY_IFDEF)"
  my @lines = @_;
  my $failed = 0;
  foreach my $line (@lines) {
    if ($line =~/^\s*#ifdef/) {
      $failed++;
    }
    elsif($line =~/^\s*#ifndef/) {
      $failed++;
    }
  }
  return $failed;
}

sub cpp_comment {
  # C pre-processor directives should not be intermingled with
  # fortran style comments
  my @lines = @_;
  my $failed = 0;
  my @comments = ();
  foreach my $line (@lines) {
    # is this an #if statement?
    if (($line =~/^\s*#if /) || ($line =~/^\s*#elif /)){
      # does this ifdef have a ! in it?
      if ($line =~/!/){
	# split the possible regions (ignoring the 0th)
	# and loop over to check each one in turn
	@comments = split/!/, $line,-1;
	splice(@comments,0,1);
	foreach my $comment (@comments) {
	  # must be a recognisable CPP directive                                                                                                                                                            
	  if ($comment !~ /(^\s*\(?\s*defined)|(^=\s*[0-9])/){
	    $failed++;
	  }
	}
      }
    }
    # is this an #else?
    elsif ($line =~/^\s*#else\s*!/){
      $failed++;
    }
    # is this an #endif?
    elsif ($line =~/^\s*#endif\s*!/){
      $failed++;
    }
    # is this an #include?
    elsif ($line =~/^\s*#include[^!]+!/){
      $failed++;
    }
  }
  return $failed;
}

sub obsolescent_fortran_intrinsic {
  my @lines = @_;

  my $failed = 0;

# Iterate over lines and keywords
  foreach my $line (@lines) {
    my @keywords_to_check = get_archaic_fortran_keywords();
# Remove comments unless they're OpenMP commands
    if ($line =~/![^\$]/) {
        $line =~s/![^\$].*//g;
    }
    $line = remove_quoted($line);

    next unless $line;
    next unless $line =~/\S/;   # If line empty, try the next
    foreach my $keyword (@keywords_to_check) {
# If the keyword is present on the line
      if ($line =~/\W$keyword\W/i or $line =~/^$keyword\W/i or
          $line =~/\%$keyword\W /i or $line =~/\($keyword\)/i) {
# If the line contains WRITE let it pass as it's probably text
        if ($line =~/WRITE.*'/ or $line =~/WRITE.*"/ ) {
          next;
        }
# Ignore CPP
	unless ($line =~/^\s*#/) {
	  $extra_error_information{$keyword}++;
	  $failed++;
        }
      }
    }
  }
  return $failed;
}
sub c_openmp_define_pair_thread_utils {
  my @lines = c_sanitise_lines(@_);

  my $failed = 0;
  for (my $i=0; $i < scalar @lines; $i++) {
    my $line = $lines[$i];

    # match ifdef and defined style for _OPENMP
    if ($line =~/^\s*#(el)?if.*defined\(_OPENMP\)/) {
      # fail if _OPENMP is not the first defined() test, or it is not
      # followed by UM_USE_C_OPENMP_VIA_THREAD_UTILS
      if ($line !~ /^\s*#(el)?if\s*!?defined\(_OPENMP\)\s*&&\s*!?defined\(UM_USE_C_OPENMP_VIA_THREAD_UTILS\)/ ) {
        $failed++;
      }
    }
  }

  return $failed;
}

sub c_openmp_define_no_combine {
  my @lines = c_sanitise_lines(@_);

  my $failed = 0;
  for (my $i=0; $i < scalar @lines; $i++) {
    my $line = $lines[$i];

    # fail if we match defined(_OPENMP) + at least two other defined()
    if ($line =~/^\s*#(el)?if\s*defined\(_OPENMP\)(.*?!?defined\(\w+\)){2,}/) {
        $failed++;
    }
  }

  return $failed;
}

sub c_openmp_define_not {
  my @lines = c_sanitise_lines(@_);

  my $failed = 0;
  for (my $i=0; $i < scalar @lines; $i++) {
    my $line = $lines[$i];

    # fail if we match !defined(_OPENMP)
    if ($line =~/^\s*#(el)?if.*!defined\(_OPENMP\)/) {
        $failed++;
    }
  }

  return $failed;
}

sub c_ifdef_defines {
  my @lines = c_sanitise_lines(@_);

  my $failed = 0;
  for (my $i=0; $i < scalar @lines; $i++) {
    my $line = $lines[$i];

    # fail if we match #ifdef or #ifndef
    if ($line =~/^\s*#if(n)?def/) {
        $failed++;
    }
  }

  return $failed;
}

sub c_protect_omp_pragma {
  my @lines = c_sanitise_lines(@_);

  # remove _OPENMP if-def protected lines.
  # As #ifs may be nested, successivly remove all the #if blocks until
  # none are remaining.
  for (my $i = scalar @lines; $i > 0; $i--) {
    my $line = $lines[$i-1];

    # as we are going from the bottom, the first #if will be an
    # innermost one.
    if ( $line =~/^\s*#if/) {

      splice @lines, $i-1, 1, '';

      my $whipe=0;

      if ($line =~/defined\(_OPENMP\)/) {
        $whipe=1;
      }

      for (my $j = $i - 1; $j < scalar @lines; $j++) {
         my $jline = $lines[$j];

         if ($whipe==1) {
           splice @lines, $j, 1, '';
         }

         if ( $jline =~/#else/ ) {
           if ($whipe==1) {
             $whipe = 0;
           }
           splice @lines, $j, 1, '';
         }

         if ( $jline =~/#elif/ ) {
           if ($whipe==1) {
             $whipe = 0;
           }
           if ($jline =~/defined\(_OPENMP\)/) {
             $whipe=1;
           }
           splice @lines, $j, 1, '';
         }

         if ( $jline =~/#endif/ ) {
           splice @lines, $j, 1, '';
           last;
         }

      }
    }
  }

  my $failed = 0;
  for (my $i=0; $i < scalar @lines; $i++) {
    my $line = $lines[$i];

    # as we have removed all lines protected by _OPENMP,
    # any remaining pragma lines are a fail.
    if ($line =~/#pragma\s+omp/) {
      $failed++;
    }

    # as are omp includes
    if ($line =~/#include\s+(<|")omp.h(>|")/) {
      $failed++;
    }
  }

  return $failed;
}

sub c_sanitise_lines {
  my @lines = @_;

  my $entire = join("",@lines);

  #remove commented sections
  $entire =~s/\/\*(.|\n)+?(\*\/)//g;

  # Sort out continuation lines
  $entire =~s/\\\s*\n//g;

  # standardise format for defined(<DEF>) style tests
  $entire =~s/defined\s*?\(?\s*?(\w+)[^\S\n]*\)?([|&><*+%^$()\/\-\s])/defined($1) $2/g;

  @lines = split/\n/,$entire;

  return @lines;
}

sub line_trail_whitespace {
  my @lines = @_;
  my $failed = 0;

  foreach my $line (@lines) {

    $line =~ s/\n//g;

    # Fail if there are whitespace characters at the end of a line.
    if ($line =~/\s+$/) {
      $failed++ ;
      $line = "\n    '$line'";
      $extra_error_information{$line}++;
    }
  }
  return $failed;
}

sub c_integral_format_specifiers {
  my @lines = @_;
  my $failed = 0;

  my @fixed_width_size = (
    '8',
    '16',
    '32',
    '64'
  );

  my @fixed_width_type = (
    'MAX',
    'PTR'
  );

  my @fixed_prefix = (
    'PRI',
    'SCN'
  );

  my @fixed_suffix = (
    '',
    'FAST',
    'LEAST'
  );

  my @print_style = (
    'd',
    'i',
    'u',
    'o',
    'x',
    'X'
  );

  # Exact numerical width style (e.g. PRIdFAST64)
  foreach my $line (@lines) {
    foreach my $fwpre (@fixed_prefix) {
      foreach my $fwps (@print_style) {
        foreach my $fwsz (@fixed_width_size) {
          foreach my $fwsfx (@fixed_suffix) {
            # Fail if format specifier immediately follows or proceeds a " character
            if ($line =~/"${fwpre}${fwps}${fwsfx}${fwsz}/) {
              $failed++ ;
              chomp($line);
              $line = "\n    '$line'";
              $extra_error_information{$line}++;
            } elsif ($line =~/${fwpre}${fwps}${fwsfx}${fwsz}"/) {
              $failed++ ;
              chomp($line);
              $line = "\n    '$line'";
              $extra_error_information{$line}++;
            }
          }
        }
      }
    }
  }

  # Style defining the width by type (e.g. SCNuMAX)
  foreach my $line (@lines) {
    foreach my $fwpre (@fixed_prefix) {
      foreach my $fwps (@print_style) {
        foreach my $fwt (@fixed_width_type) {
          # Fail if format specifier immediately follows or proceeds a " character
          if ($line =~/"${fwpre}${fwps}${fwt}/) {
            $failed++ ;
            chomp($line);
            $line = "\n    '$line'";
            $extra_error_information{$line}++;
          } elsif ($line =~/${fwpre}${fwps}${fwt}"/) {
            $failed++ ;
            chomp($line);
            $line = "\n    '$line'";
            $extra_error_information{$line}++;
          }
        }
      }
    }
  }

  return $failed;
}

sub c_final_newline {
  my @lines = @_;
  my $failed = 0;

  my $line = $lines[-1];
  my $fchar = substr $line, -1;

  # Fail if the final line does not end with a newline character
  if (ord $fchar != 10) {
    $failed++ ;
  }

  return $failed;
}

1;
