/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   ft_memmove.c                                       :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: igomez <igomez@student.42.fr>              +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2014/11/05 12:08:53 by igomez            #+#    #+#             */
/*   Updated: 2014/11/12 17:12:51 by igomez           ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

#include "libft.h"

void	*ft_memmove(void *dst, const void *src, size_t len)
{
	size_t		i;

	if (!src || !dst)
		return (NULL);
	if (src < dst)
	{
		i = 0;
		while (i < len)
		{
			((char *)dst)[len - 1 - i] = ((char *)src)[len - 1 - i];
			i++;
		}
	}
	else
		ft_memcpy(dst, src, len);
	return (dst);
}
