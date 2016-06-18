/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   ft_str_cleanpath.c                                 :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: igomez <igomez@student.42.fr>              +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2015/01/20 18:41:10 by igomez            #+#    #+#             */
/*   Updated: 2015/01/20 18:41:11 by igomez           ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

#include "libft.h"

char	*ft_clean_double_slash(char *s)
{
	size_t	i;
	size_t	j;

	if (!s)
		return (s);
	i = 0;
	j = 0;
	while ((s[i] = s[j]))
	{
		if (!(s[j] == '/' && s[j + 1] == '/'))
		{
			i += 1;
			j += 1;
		}
		else
			j += 1;
	}
	if ((i > 1) && ft_strequ(s + i - 1, "/"))
		s[i - 1] = 0;
	return (s);
}

char	*ft_clean_same_dir(char *s)
{
	size_t	i;
	size_t	j;

	if (!s)
		return (s);
	i = 0;
	j = 0;
	while ((s[i] = s[j]))
	{
		if (!ft_strnequ(s + j, "/./", 3))
		{
			i += 1;
			j += 1;
		}
		else
			j += 2;
	}
	if ((i > 2) && ft_strequ(s + i - 2, "/."))
		s[i - 2] = 0;
	return (s);
}

char	*ft_clean_parent(char *s)
{
	char	*cursor;
	size_t	i;
	size_t	j;

	if (!s || *s != '/')
		return (s);
	i = 0;
	while (s[i])
	{
		if (ft_strnequ(s + i, "/../", 4))
		{
			s[i] = 0;
			if (!(cursor = ft_strrchr(s, '/')))
				cursor = s;
			j = 0;
			while ((cursor[j] = s[i + j + 3]))
				j++;
			i = cursor - s;
		}
		else
			i++;
	}
	return (s);
}

char	*ft_clean_trailing_dotdot(char *s)
{
	char	*cursor;

	if (!s || *s != '/')
		return (s);
	if (ft_strequ(s, "/.."))
		s[1] = 0;
	else if ((cursor = ft_strrchr(s, '/')) && ft_strequ(cursor, "/.."))
	{
		*cursor = 0;
		cursor = ft_strrchr(s, '/');
		cursor[cursor == s ? 1 : 0] = 0;
	}
	return (s);
}

char	*ft_str_cleanpath(char *s)
{
	if (!s)
		return (s);
	ft_clean_double_slash(s);
	ft_clean_same_dir(s);
	if (*s == '/')
	{
		ft_clean_parent(s);
		ft_clean_trailing_dotdot(s);
	}
	return (s);
}
