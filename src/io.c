
// Prints program logo
void logo(int version)
{
	border_print();
	printf(
	"                   __   __ ___________                 _                        \n"
	"                   \\ \\ / //  ___| ___ \\               | |                       \n"
	"                    \\ V / \\ `--.| |_/ / ___ _ __   ___| |__                     \n"
	"                    /   \\  `--. \\ ___ \\/ _ \\ '_ \\ / __| '_ \\                    \n"
	"                   / /^\\ \\/\\__/ / |_/ /  __/ | | | (__| | | |                   \n"
	"                   \\/   \\/\\____/\\____/ \\___|_| |_|\\___|_| |_|                   \n\n"
    );
	border_print();
	center_print("Developed at Argonne National Laboratory", 79);
	char v[100];
	sprintf(v, "Version: %d", version);
	center_print(v, 79);
	border_print();
}

// Prints Section titles in center of 80 char terminal
void center_print(const char *s, int width)
{
	int length = strlen(s);
	int i;
	for (i=0; i<=(width-length)/2; i++) {
		fputs(" ", stdout);
	}
	fputs(s, stdout);
	fputs("\n", stdout);
}

void border_print(void)
{
	printf(
	"==================================================================="
	"=============\n");
}

// Prints comma separated integers - for ease of reading
void fancy_int( int a )
{
    if( a < 1000 )
        printf("%d\n",a);

    else if( a >= 1000 && a < 1000000 )
        printf("%d,%03d\n", a / 1000, a % 1000);

    else if( a >= 1000000 && a < 1000000000 )
        printf("%d,%03d,%03d\n", a / 1000000, (a % 1000000) / 1000, a % 1000 );

    else if( a >= 1000000000 )
        printf("%d,%03d,%03d,%03d\n",
               a / 1000000000,
               (a % 1000000000) / 1000000,
               (a % 1000000) / 1000,
               a % 1000 );
    else
        printf("%d\n",a);
}
