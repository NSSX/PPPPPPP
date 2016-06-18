/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   ft_strequ_nc.c                                     :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: igomez <igomez@student.42.fr>              +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2015/01/03 14:40:50 by igomez            #+#    #+#             */
/*   Updated: 2015/01/03 14:41:01 by igomez           ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

#include "libft.h"

int		ft_strequ_nc(char const *s1, char const *s2)
{
	int		i;

	if (!s1 || !s2)
		return (0);
	i = 0;
	while (s1[i] && ft_tolower(s1[i]) == ft_tolower(s2[i]))
		i++;
	return ((s1[i] == s2[i]) ? 1 : 0);
}
