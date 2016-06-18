/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   ft_tabdup.c                                        :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: igomez <igomez@student.42.fr>              +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2015/01/09 16:44:51 by igomez            #+#    #+#             */
/*   Updated: 2015/01/09 16:44:52 by igomez           ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

#include "libft.h"

char	**ft_tabdup(char **src)
{
	size_t	i;
	char	**dst;

	if (!src || !(dst = ft_tabnew(ft_tablen(src))))
		return (NULL);
	i = 0;
	while (src[i])
	{
		if (!(dst[i] = ft_strdup(src[i])))
		{
			ft_tabdel(&dst);
			return (NULL);
		}
		i++;
	}
	return (dst);
}
