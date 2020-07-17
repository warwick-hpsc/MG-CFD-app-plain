#include <stdio.h>

#include "progress.h"

void show_progress_bar() {
	int width=40;
	printf("[");
	for (int i=0; i<width; i++) {
		printf(" ");
	}
	printf("] 0%%\r");
	fflush(stdout);
}

void update_progress_bar(float percent) {
	if (percent > 1.0) percent = 1.0;
	if (percent < 0.0) percent = 0.0;
	int width=40;
	printf("[");
	int pos = width*percent;
	for (int i=0; i<width; i++) {
		if (i < pos) printf("=");
		else if (i == pos) printf(">");
		else printf(" ");
	}
	printf("] %d%%\r", (int)(percent*100));
	fflush(stdout);
}

void stop_progress_bar() {
	printf("\n");
}
