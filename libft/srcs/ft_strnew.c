/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   ft_strnew.c                                        :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: igomez <igomez@student.42.fr>              +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2014/11/06 13:29:44 by igomez            #+#    #+#             */
/*   Updated: 2014/11/06 15:13:55 by igomez           ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

#include "libft.h"
#include <stdlib.h>

char	*ft_strnew(size_t size)
{
	char	*dst;

	if (!(dst = (char *)malloc(sizeof(char) * (size + 1))))
		return (NULL);
	ft_memset(dst, '\0', size + 1);
	return (dst);
}
