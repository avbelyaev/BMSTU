#include <stdio.h>

int count;

int check_for_power_of_2 (int _value){
        int v1 = _value;
	if (v1 == 0) return 0;//0 точно не степень двойки
	//циклически сдвигаем число вправо, пока не останется самая старшая единица

	while(v1 != 1) { //пока число больше единицы
		if ((v1 & 1) != 0) return 0;//нашлась лишняя единица
		v1 = v1 >> 1;//сдвиг
	}
	return 1;
}

int main()
{
        int N, i;
	scanf ("%d", &N);
	int array[N];

        for (i = 0; i < N; i++) scanf("%d", &array[i]);

	int mask_value;
	//делаем маску из единиц
	unsigned int mask_limit = 0;
	for (i = 0; i != N; i++)
		mask_limit = mask_limit | (1 << i);

        count = 0;
	//цикл от 0 до максимального значения числа
	for (mask_value = 0; mask_value <= mask_limit; mask_value++){
		int sum = 0;
		int k = 0;
		//каждое значение mask_value уникально
		//след-но уникальна и последовательность нулей и единиц (из которых онро состоит)
		for (k = 0; k != N; k++)
			if (((mask_value >> k) & 1) == 1){//берем только единицы
				//берем числа (по индексу) из исходного массива
				//суммируем
				sum += array[k];
				//printf("%d+",array[k]);
			}
		//проверяем сумму на степень двойки
		if (check_for_power_of_2(sum) == 1) {
			//printf(" <power of two");
			count++;
		}
		//printf("\n");
	}
        printf("%d", count);
	return 0;
}