void revarray(void* base, unsigned long len, unsigned long width)
{
	char* left = (char*)base; /*начало массива*/
	char* right = left + ((len-1)*width); /*конец массива*/
	char temp_array[255];

	while(left + width <= right){
		memcpy(temp_array,left,width);/*копирует  width байт из left в temp_arr*/
		memcpy(left,right,width);
		memcpy(right,temp_array,width);
		left += width;
		right -= width;
	}


}