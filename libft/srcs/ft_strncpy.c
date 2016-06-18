/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   ft_strncpy.c                                       :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: igomez <igomez@student.42.fr>              +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2014/11/04 15:21:35 by igomez            #+#    #+#             */
/*   Updated: 2014/11/06 15:44:56 by igomez           ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

#include "libft.h"

char *ft_strncpy(char *dst, const char *src, size_t n)
{
	size_t	i;
	size_t	src_len;

	if (!dst || !src)
		return (NULL);
	src_len = ft_strlen(src);
	i = 0;
	while (i < n)
	{
		dst[i] = (i < src_len) ? src[i] : '\0';
		i++;
	}
	return (dst);
}
