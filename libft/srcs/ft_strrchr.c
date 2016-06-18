/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   ft_strrchr.c                                       :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: igomez <igomez@student.42.fr>              +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2014/11/04 18:19:39 by igomez            #+#    #+#             */
/*   Updated: 2014/11/07 20:10:25 by igomez           ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

#include "libft.h"

char *ft_strrchr(const char *s, int c)
{
	size_t	i;
	size_t	str_len;

	if (!s)
		return (NULL);
	str_len = ft_strlen(s);
	if ((char)c == 0)
		return ((char *)(s + str_len));
	i = 0;
	while (i < str_len)
	{
		if (s[str_len - i - 1] == (char)c)
			return ((char *)(s + str_len - i - 1));
		i++;
	}
	return (NULL);
}
