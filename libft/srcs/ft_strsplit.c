/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   ft_strsplit.c                                      :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: igomez <igomez@student.42.fr>              +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2014/11/07 10:47:53 by igomez            #+#    #+#             */
/*   Updated: 2014/11/08 15:33:24 by igomez           ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

#include "libft.h"
#include <stdlib.h>

static size_t		ft_offset_nextblock(char const *s, char c)
{
	size_t		j;

	j = 0;
	while (s[j] && ((s[0] == c && s[j] == c) || (s[0] != c && s[j] != c)))
		j++;
	return (j);
}

static size_t		ft_getnbword(char const *s, char c)
{
	size_t		i;
	size_t		j;
	size_t		n;

	i = 0;
	n = 0;
	while (s[i])
	{
		j = ft_offset_nextblock(s + i, c);
		if (s[i] != c)
			n++;
		i = i + j;
	}
	return (n);
}

char				**ft_strsplit(char const *s, char c)
{
	char	**res;
	size_t	i;
	size_t	j;
	size_t	n;

	if (!s || !(res = ft_tabnew(ft_getnbword(s, c))))
		return (NULL);
	i = 0;
	n = 0;
	while (s[i])
	{
		j = ft_offset_nextblock(s + i, c);
		if (s[i] != c)
		{
			if (!(res[n++] = ft_strsub(s, i, j)))
			{
				ft_tabdel(&res);
				return (NULL);
			}
		}
		i = i + j;
	}
	return (res);
}
