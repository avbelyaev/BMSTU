//ИУ4_какое-то_ДЗ_сем2

#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <time.h>
#define MAXLEN 40
#define READLEN 300

//changelog:
            //minor bug fixes
            //memory allocated as much as prog needs (prev size:1000000)
            //added time stamps
            //changed bubblesort -> shellsort

typedef struct student {
    char first[MAXLEN];
    char second[MAXLEN];
    char univer[MAXLEN];
    char fac[MAXLEN];
    char dep[MAXLEN];
    int age, course, group;
    char city[MAXLEN];
} stud_tag;

struct town {
    char name[MAXLEN];
    int stud, bach, mag;
};

void swap (stud_tag *i, stud_tag *j) {
    stud_tag temp = *i;
    *i = *j;
    *j = temp;
}

int main()
{
    FILE *fp = fopen("students_1.csv", "r");
    FILE *ff = fopen("result_1.txt", "w");

    if (NULL == fp) {
        printf("cannot open file\n");
        exit(1);
    }

    int i = 0, j = 0, stud_num = 0, memory_size = 0, bachelors = 0, magisters = 0, eldest, city_num = 0, mark, k = 0, gap;
    double start_time, c_t1, c_t2, c_t3, c_t4, c_t5, finish_time;
    char str[READLEN];

    while (!feof(fp)) if (fgets(str, READLEN, fp)) memory_size++;  //чтобы не выделять лишней памяти считаем сначала сколько всего строк в файле
    memory_size++;
    rewind(fp);

    struct student* man;
    man = (struct student *)malloc(sizeof(struct student)*memory_size);
    printf("memory allocated\n");

    start_time = clock();
    printf("\tstart:%.2lfms\n", start_time);  //здесь именно текущее время

    printf("reading file: ");
    fscanf (fp, "%*[a-zA-Z];%*[a-zA-Z];%*[a-zA-Z];%*[a-zA-Z];%*[a-zA-Z];%*[a-zA-Z];%*[a-zA-Z];%*[a-zA-Z];%*[a-zA-Z];");
    while (9 == fscanf(fp, "%[^;];%[^;];%[^;];%[^;];%[^;];%d;%d;%d;%[^;];\n", man[i].first, man[i].second, man[i].univer,
                man[i].fac, man[i].dep, &man[i].age, &man[i].course, &man[i].group, man[i].city)) {
        if (5 > man[i].course) {
            bachelors++;
        } else {
            magisters++;
        }
        i++;
    }
    fclose(fp);
    printf("done\n");
    c_t1 = clock();
    printf("\treading duration:%.2lfms\n", c_t1 - start_time);  //а здесь именно продолжительность процесса

    //===========LEVEL1===========

    stud_num = i;
    fprintf(ff, "Total number of students:%d\nBachelors:%d / Magisters:%d\n", stud_num, bachelors, magisters);

    eldest = man[0].age;
    for (i = 0; i < stud_num; i++) if (man[i].age > eldest) eldest = man[i].age;
    //Bubblesort
    /*for (i = 0; i < stud_num - 1; i++) {
        for (j = i + 1; j < stud_num; j++) {
            if (0 == strcmp(man[i].second, man[j].second) && 0 < strcmp(man[i].first, man[j].first)) {
                swap(&man[i], &man[j]);
            } else {
                if (0 < strcmp(man[i].second, man[j].second)) swap(&man[i], &man[j]);
            }
        }
    }*/
    //sHELLsort
    for (gap = stud_num/2; gap > 0; gap /= 2) //wiki
        for (i = gap; i < stud_num; i++) 
		for (j = i; j >= gap && ((0 == strcmp(man[j-gap].second, man[j].second) && 0 < strcmp(man[j-gap].first, man[j].first)) || (0 < strcmp(man[j-gap].second, man[j].second))); j -= gap) swap (&man[j], &man[j - gap]);

    fprintf(ff, "\nAlphabetical sorting:\n");
    for (i = 0; i < stud_num; i++) fprintf(ff, "\t%s %s\n", man[i].second, man[i].first);
    printf(">>>Level 1 Finished\n");
    printf("sorting done\n");
    c_t2 = clock();
    printf("\tsorting duration:%.2lfms\n", c_t2 - c_t1);

    //===========LEVEL2===========

    fprintf(ff, "\nEldest student's age:%d\n\nEldest students:\n", eldest);
    for (j = 1, i = 0; i < stud_num; i++) if (man[i].age == eldest) fprintf(ff, "%d)%s %s, %d\n\t%s %s\n\tfaculty:%s course:%d group:%d\n\t%s\n\n",
                                                                            j++ , man[i].second, man[i].first, man[i].age, man[i].univer, man[i].fac,
                                                                            man[i].dep, man[i].course, man[i].group, man[i].city);
    printf(">>>Level 2 Finished\n");
    c_t3 = clock();
    printf("\t'eldest' duration:%.2lfms\n", c_t3 - c_t2);

    //===========LEVEL3===========

    struct town* cities;
    cities = (struct town *)malloc(sizeof(struct town)*memory_size);
    printf("memory2 allocated\n");

    for (i = 0; i < stud_num; i++) {
        for (mark = 0, j = 0; j < i; j++) {
            if (0 == strcmp(man[i].city, man[j].city)) {
                mark++;
                break;
            }
        }
        if (0 == mark) {
            strcpy(cities[k].name, man[i].city);
            city_num++;
            k++;
        }
    }
    printf("array filled\n");
    fprintf(ff, "\n\nNumber of cities:%d\n", city_num);
    c_t4 = clock();
    printf("\tfill duration:%.2lfms\n", c_t4 - c_t3);

    for (i = 0; i < city_num; i++) { cities[i].stud = 0; cities[i].bach = 0; cities[i].mag = 0; }

    for (i = 0; i < city_num; i++) {
        for (j = 0; j < stud_num; j++) {
            if (0 == strcmp(cities[i].name, man[j].city)) {
                if (5 > man[j].course) {
                    cities[i].bach++;
                } else {
                    cities[i].mag++;
                }
                cities[i].stud++;
            }
        }
    }
    printf("courses distributed\n");
    c_t5 = clock();
    printf("\tdistr duration:%.2lfms\n", c_t5 - c_t4);
    fprintf(ff, "Distribution by cities:\n");
    for (i = 0; i < city_num; i++) fprintf(ff, "%s:\n\ttotal students:%d\n\tbachelors:%d / magisters:%d\n\n", cities[i].name, cities[i].stud, cities[i].bach, cities[i].mag);

    int most_stud = cities[0].stud, most_bachs = cities[0].bach, most_mags = cities[0].mag;

    for (i = 1; i < city_num; i++) {
        if (cities[i].stud > most_stud) most_stud = cities[i].stud;
        if (cities[i].bach > most_bachs) most_bachs = cities[i].bach;
        if (cities[i].mag > most_mags) most_mags = cities[i].mag;
    }

    printf("most_courses counted\n");
    for (i = 0; i < city_num; i++) {
        if (cities[i].stud == most_stud) fprintf(ff, "Most students from %s\n", cities[i].name);
        if (cities[i].bach == most_bachs) fprintf(ff, "Most bachelors from %s\n", cities[i].name);
        if (cities[i].mag == most_mags) fprintf(ff, "Most magisters from %s\n", cities[i].name);
    }

    fprintf(ff, "\nAverage num of students:%d\nAvg num of bachelors:%d / magisters:%d\n\n", stud_num / city_num, bachelors / city_num, magisters / city_num);
    printf(">>>Level 3 Finished\n");
    finish_time = clock();
    printf("\tTotal duration:%.2lfms\n", finish_time - start_time);

    fclose(ff);
    return 0;
//A.B
}