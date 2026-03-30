#include "utility.h"
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#define SUMMARY_FIELDS 5

typedef struct BlockSummary {
  unsigned long long pair_count;
  unsigned long long first_odd;
  unsigned long long last_odd;
  int first_is_prime;
  int last_is_prime;
} BlockSummary;

static unsigned long long last_odd_in_range(unsigned long long range_end) {
  if (range_end < 3ULL) {
    return 0ULL;
  }

  return (range_end & 1ULL) ? range_end : (range_end - 1ULL);
}

static int is_prime(unsigned long long value) {
  unsigned long long divisor = 0ULL;

  if (value < 2ULL) {
    return 0;
  }

  if (value == 2ULL) {
    return 1;
  }

  if (value % 2ULL == 0ULL) {
    return 0;
  }

  for (divisor = 3ULL; divisor <= value / divisor; divisor += 2ULL) {
    if (value % divisor == 0ULL) {
      return 0;
    }
  }

  return 1;
}

static BlockSummary process_block(unsigned long long first_odd,
                                  unsigned long long odd_count) {
  BlockSummary summary = {0, 0, 0, 0, 0};
  unsigned long long current = 0;
  unsigned long long index = 0;
  int previous_is_prime = 0;

  if (odd_count == 0ULL) {
    return summary;
  }

  summary.first_odd = first_odd;
  summary.last_odd = first_odd + 2ULL * (odd_count - 1ULL);

  current = first_odd;
  for (index = 0ULL; index < odd_count; ++index, current += 2ULL) {
    int current_is_prime = is_prime(current);

    if (index == 0ULL) {
      summary.first_is_prime = current_is_prime;
    } else if (previous_is_prime && current_is_prime) {
      summary.pair_count++;
    }

    previous_is_prime = current_is_prime;
  }

  summary.last_is_prime = previous_is_prime;
  return summary;
}

int main(int argc, char **argv) {
  Args ins__args;
  struct timeval ins__tstart, ins__tstop;
  const unsigned long long range_start = 2ULL;
  unsigned long long range_end = 0ULL;
  unsigned long long first_odd = 0ULL;
  unsigned long long last_odd = 0ULL;
  unsigned long long total_odds = 0ULL;
  unsigned long long base_odds = 0ULL;
  unsigned long long extra_odds = 0ULL;
  unsigned long long local_count = 0ULL;
  unsigned long long local_offset = 0ULL;
  unsigned long long local_start = 0ULL;
  unsigned long long total_pairs = 0ULL;
  unsigned long long local_summary[SUMMARY_FIELDS] = {0ULL, 0ULL, 0ULL, 0ULL,
                                                      0ULL};
  int myrank = 0;
  int nproc = 0;

  parseArgs(&ins__args, &argc, argv);
  range_end = (ins__args.arg >= 2L) ? (unsigned long long)ins__args.arg : 1ULL;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  last_odd = last_odd_in_range(range_end);
  if (last_odd != 0ULL) {
    first_odd = 3ULL;
    total_odds = ((last_odd - first_odd) / 2ULL) + 1ULL;
  }

  base_odds = total_odds / (unsigned long long)nproc;
  extra_odds = total_odds % (unsigned long long)nproc;
  local_count = base_odds + (((unsigned long long)myrank < extra_odds) ? 1ULL : 0ULL);
  local_offset = base_odds * (unsigned long long)myrank +
                 (((unsigned long long)myrank < extra_odds) ? (unsigned long long)myrank
                                                            : extra_odds);

  MPI_Barrier(MPI_COMM_WORLD);
  if (myrank == 0) {
    gettimeofday(&ins__tstart, NULL);
  }

  if (local_count > 0ULL) {
    BlockSummary summary;

    local_start = first_odd + 2ULL * local_offset;
    summary = process_block(local_start, local_count);

    local_summary[0] = summary.pair_count;
    local_summary[1] = summary.first_odd;
    local_summary[2] = summary.last_odd;
    local_summary[3] = (unsigned long long)summary.first_is_prime;
    local_summary[4] = (unsigned long long)summary.last_is_prime;
  }

  if (myrank == 0) {
    unsigned long long previous_last_odd = 0ULL;
    int previous_last_is_prime = 0;
    int rank = 0;

    total_pairs += local_summary[0];
    if (local_summary[1] != 0ULL) {
      previous_last_odd = local_summary[2];
      previous_last_is_prime = (int)local_summary[4];
    }

    for (rank = 1; rank < nproc; ++rank) {
      unsigned long long remote_summary[SUMMARY_FIELDS];

      MPI_Recv(remote_summary, SUMMARY_FIELDS, MPI_UNSIGNED_LONG_LONG, rank, 0,
               MPI_COMM_WORLD, MPI_STATUS_IGNORE);

      total_pairs += remote_summary[0];

      if (previous_last_odd != 0ULL && remote_summary[1] != 0ULL &&
          previous_last_odd + 2ULL == remote_summary[1] &&
          previous_last_is_prime && remote_summary[3] != 0ULL) {
        total_pairs++;
      }

      if (remote_summary[1] != 0ULL) {
        previous_last_odd = remote_summary[2];
        previous_last_is_prime = (int)remote_summary[4];
      }
    }

    gettimeofday(&ins__tstop, NULL);
    printf("Twin primes in [%llu, %llu]: %llu\n", range_start, range_end,
           total_pairs);
    ins__printtime(&ins__tstart, &ins__tstop, ins__args.marker);
  } else {
    MPI_Send(local_summary, SUMMARY_FIELDS, MPI_UNSIGNED_LONG_LONG, 0, 0,
             MPI_COMM_WORLD);
  }

  MPI_Finalize();
  return 0;
}
