#include "Utility.h"

double randomPercentage(){
    return (double)rand()/(double)RAND_MAX;
}

void appendToStr(char **str, char *toAppend, int *count, int *capacity){
    if((*count) == (*capacity)){
        // resize
        *capacity *= 2;
        char *tmp = (char *)realloc((*str), sizeof(char) * (*capacity));
        if(!tmp){
            printf("Failed to reallocate memory for matrix string!\n");
            exit(-1);
        }
        *str = tmp;
    }
    
    int i;
    for(i = 0; i < strlen(toAppend); i++){
        (*str)[(*count)++] = toAppend[i];
        if((*count) == (*capacity)){
            // resize
            *capacity *= 2;
            char *tmp = (char *)realloc((*str), sizeof(char) * (*capacity));
            if(!tmp){
                printf("Failed to reallocate memory for matrix string!\n");
                exit(-1);
            }
            *str = tmp;
        }
        (*str)[*count] = '\0';
    }
}

char * readLine(FILE **file){
    char c;
    int endOfLine = 0;
    int lineCapactiy = 16;
    char *line = (char*)malloc(sizeof(char) * lineCapactiy);
    char *tmp;
    line[0] = '\0';

    c = getc(*file);
    while(c != '\n' && c != EOF){
        line[endOfLine++] = c;
        line[endOfLine] = '\0';
        if(endOfLine == lineCapactiy - 1){
            lineCapactiy *= 2;
            tmp = (char *)realloc(line, sizeof(char) * lineCapactiy);
            if(tmp == NULL){
                printf("Failed to reallocate space for line.");
                free(tmp);
                exit(1);
            }
            line = tmp;
        }
        c = getc(*file);
    }

    return line;
}
