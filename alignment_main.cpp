// crt_aligned_malloc.c

#include <malloc.h>
#include <stdio.h>
#include <iostream>

int main() {
	void    *ptr;
	size_t  alignment,
		off_set;

	// Note alignment should be 2^N where N is any positive int.
	alignment = 16;
	off_set = 5;

	// Using _aligned_malloc
	ptr = _aligned_malloc(100, alignment);
	if (ptr == NULL)
	{
		printf_s("Error allocation aligned memory.");
		return -1;
	}
	if (((unsigned long long)ptr % alignment) == 0)
		printf_s("This pointer, %p, is aligned on %Iu\n",
		ptr, alignment);
	else
		printf_s("This pointer, %p, is not aligned on %Iu\n",
		ptr, alignment);

	// Using _aligned_realloc
	ptr = _aligned_realloc(ptr, 200, alignment);
	if (((unsigned long long)ptr % alignment) == 0)
		printf_s("This pointer, %p, is aligned on %Iu\n",
		ptr, alignment);
	else
		printf_s("This pointer, %p, is not aligned on %Iu\n",
		ptr, alignment);
	_aligned_free(ptr);

	// Using _aligned_offset_malloc
	ptr = _aligned_offset_malloc(200, alignment, off_set);
	if (ptr == NULL)
	{
		printf_s("Error allocation aligned offset memory.");
		return -1;
	}
	if (((((unsigned long long)ptr) + off_set) % alignment) == 0)
		printf_s("This pointer, %p, is offset by %Iu on alignment of %Iu\n",
		ptr, off_set, alignment);
	else
		printf_s("This pointer, %p, does not satisfy offset %Iu "
		"and alignment %Iu\n", ptr, off_set, alignment);

	// Using _aligned_offset_realloc
	ptr = _aligned_offset_realloc(ptr, 200, alignment, off_set);
	if (ptr == NULL)
	{
		printf_s("Error reallocation aligned offset memory.");
		return -1;
	}
	if (((((unsigned long long)ptr) + off_set) % alignment) == 0)
		printf_s("This pointer, %p, is offset by %Iu on alignment of %Iu\n",
		ptr, off_set, alignment);
	else
		printf_s("This pointer, %p, does not satisfy offset %Iu and "
		"alignment %Iu\n", ptr, off_set, alignment);

	// Note that _aligned_free works for both _aligned_malloc
	// and _aligned_offset_malloc. Using free is illegal.
	_aligned_free(ptr);



	std::system("pause");
}