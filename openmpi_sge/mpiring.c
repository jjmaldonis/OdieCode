int main(int argc, char *argv[])

{

	MPI_Comm myComm; /* intra-communicator of local sub-group */

	MPI_Comm myFirstComm; /* inter-communicators */

	MPI_Comm mySecondComm;

	int membershipKey;

	int rank;

	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	...

		/* User code must generate membershipKey in the range [0, 1, 2] */

		membershipKey = rank % 3;

	/* Build intra-communicator for local sub-group */

	MPI_Comm_split(MPI_COMM_WORLD, membershipKey, rank, &myComm);

	/* Build inter-communicators. Tags are hard-coded. */

	if (membershipKey == 0)

	{ /* Group 0 communicates with groups 1 and 2. */

		MPI_Intercomm_create( myComm, 0, MPI_COMM_WORLD, 1,

				1, &myFirstComm);

		MPI_Intercomm_create( myComm, 0, MPI_COMM_WORLD, 2,

				2, &mySecondComm);

	}

	else if (membershipKey == 1)

	{ /* Group 1 communicates with groups 0 and 2. */

		MPI_Intercomm_create( myComm, 0, MPI_COMM_WORLD, 0,

				1, &myFirstComm);

		MPI_Intercomm_create( myComm, 0, MPI_COMM_WORLD, 2,

				12, &mySecondComm);

	}

	else if (membershipKey == 2)

	{ /* Group 2 communicates with groups 0 and 1. */

		MPI_Intercomm_create( myComm, 0, MPI_COMM_WORLD, 0,

				2, &myFirstComm);

		MPI_Intercomm_create( myComm, 0, MPI_COMM_WORLD, 1,

				12, &mySecondComm);

	}

	/* Do some work ... */

	/* Then free communicators before terminating... */

	MPI_Comm_free(&myFirstComm);

	MPI_Comm_free(&mySecondComm);

	MPI_Comm_free(&myComm);

	MPI_Finalize();

	return 0;

}
